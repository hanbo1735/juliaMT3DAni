"""
 module `MT1DFwdAni` defines routines to solve MT 1D anisotropic forward
 modeling problem.

 - by Bo Han, Oct 2016

"""

module MT1DFwdAni

export mt1DAniAnalyticField, mt1DAniAnalyticImp


"""
`mt1DAniAnalyticField` computes analytic fields for 1D general anisotropic model.
 The solution is basd on field propagation formulas (see the document by Yuguo Li).
 This function is basically wrapped from Yuguo's FORTRAN counterpart.

Input:
    freq    :: Float64  - a single frequency value.
    sigma   :: Array    - nLayer*6, conductivity tensor of layered model.
    zNode   :: Vector   - nLayer+1, depth of top of each layer.
    eTop    :: Vector   - [Ex0; Ey0], the given top boundary E-field value.
    compH   :: Bool     - whether compute the H-field or not.
    timeFac :: String   - which time factor will be used: e^{-iwt} ("neg",
                          by default) or e^{iwt} ("pos")

Output:
    Ex, Ey, Ez, Hx, Hy :: Vector{Complex}
                          - EM fields for each layer (Ex, Ey, Hx, Hy at layer
                            interfaces, Ez at layer centers).

"""
function mt1DAniAnalyticField{T<:Float64}(freq::T, sigma::Array{T}, zNode::Vector{T},
                                          eTop::Vector{Complex{T}}; compH::Bool=false,
                                          timeFac::AbstractString="neg")

    # check if layer conductivity and depth array have the correct size
    if size(sigma,1) != length(zNode)-1
        error("The dimension of layer conductivity and depth mismatch!")
    end

    # physical constant
    MU0   = 4 * pi * 1e-7
    omega = 2 * pi * freq
    iom   = 1im * omega * MU0    # Time dependence: e^{-iwt}

    # assume halfspace at the bottom of model
    sig    = vcat(sigma, sigma[end:end, :])
    nLayer = length(zNode)
    zLen   = diff(zNode)

    # effective azimuthal anisotropic conductivity
    sigEff = getEffSigma(sig)

    sh = (x::Complex) -> 0.5*(exp(x) - exp(-x))   # sinh()
    ch = (x::Complex) -> 0.5*(exp(x) + exp(-x))   # cosh()

    Sprod = zeros(Complex128, 4, 4, nLayer)
    Sprod[:, :, end] = eye(4)

    M = zeros(Complex{T}, 4, 4)

    xip = zeros(Complex128, nLayer)
    xim = copy(xip)

    QpNeg    = zeros(Float64, nLayer)
    QmNegRec = copy(QpNeg)

    # loop over layers to get the accumulative product of S (i.e. Sprod).
    for j = nLayer:-1:1

        Axx = sigEff[j, 1]
        Ayy = sigEff[j, 2]
        Axy = sigEff[j, 3]

        ada = Axx + Ayy
        add = Axx - Ayy

        dA12 = sqrt( add^2 + 4*Axy^2 )
        if add<0; dA12 = -dA12; end

        A1 = 0.5 * (ada + dA12)
        A2 = 0.5 * (ada - dA12)

        kp = sqrt(-iom) * sqrt(A1)
        km = sqrt(-iom) * sqrt(A2)

        # xip, xim, Qp, and Qm for each layer are saved for later use.
        xip[j] = -kp/iom
        xim[j] = -km/iom

        QpNeg[j]    = Axy == 0 ? 0 : 2*Axy/(add+dA12)         # -Qp
        QmNegRec[j] = Axy == 0 ? 0 : 0.5*(add-dA12)/Axy       # -1/Qm


        QpDerQm = QpNeg[j] * QmNegRec[j]                # Qp/Qm
        dq = 1-QpDerQm                                  # 1-Qp/Qm

        # M of the basement
        if j==nLayer

            # 1st and 3rd columns won't be used, so set them to zeros
            M = [ 0                  1     0          QmNegRec[j];
                  0           QpNeg[j]     0                    1;
                  0   -xip[j]*QpNeg[j]     0              -xim[j];
                  0             xip[j]     0   xim[j]*QmNegRec[j] ]

            continue
        end

        sp = sh(kp*zLen[j])
        sm = sh(km*zLen[j])
        cp = ch(kp*zLen[j])
        cm = ch(km*zLen[j])

        S = zeros(Complex{T}, 4, 4)
        S[1, 1] =  cp-QpDerQm*cm
        S[1, 2] = -QmNegRec[j]*(cp-cm)
        S[1, 3] =  QmNegRec[j]*(sp/xip[j]-sm/xim[j])
        S[1, 4] =  sp/xip[j]-QpDerQm*sm/xim[j]
        S[2, 1] =  QpNeg[j]*(cp-cm)
        S[2, 2] = -QpDerQm*cp+cm
        S[2, 3] =  QpDerQm*sp/xip[j]-sm/xim[j]
        S[2, 4] =  QpNeg[j]*(sp/xip[j]-sm/xim[j])
        S[3, 1] = -QpNeg[j]*(xip[j]*sp-xim[j]*sm)
        S[3, 2] =  QpDerQm*xip[j]*sp-xim[j]*sm
        S[3, 3] = -QpDerQm*cp+cm
        S[3, 4] = -QpNeg[j]*(cp-cm)
        S[4, 1] =  xip[j]*sp-QpDerQm*xim[j]*sm
        S[4, 2] = -QmNegRec[j]*(xip[j]*sp-xim[j]*sm)
        S[4, 3] =  QmNegRec[j]*(cp-cm)
        S[4, 4] =  cp-QpDerQm*cm

        S = S/dq

        Sprod[:,:,j] = S * Sprod[:,:,j+1]

    end  # nLayer


    SM = Sprod[:,:,1] * M

    # assume Ex0 & Ey0 are known, compute Cm and Dm of the bottom layer
    Ex0 = eTop[1]
    Ey0 = eTop[2]
    G = [SM[1,2]  SM[1,4]; SM[2,2]  SM[2,4]]

    CD = G \ [Ex0; Ey0]
    # CD = inv(G) * [Ex0; Ey0]
    #CD    = zeros(Complex128, 2, 1)
    #CD[1] = ( G[2,2]*Ex0 - G[1,2]*Ey0)/det(G)
    #CD[2] = (-G[2,1]*Ex0 + G[1,1]*Ey0)/det(G)

#     # Alternatively, one can assume Hx0 & Hy0 are known
#     Hx0 = 1.
#     Hy0 = 1.
#     G = [SM[3,2]  SM[3,4]; SM[4,2]  SM[4,4]]
#     CD = inv(G) * [Hx0; Hy0]

    # the full C of the bottom layer
    Cbot = [0; CD[1]; 0; CD[2]]

    MC = M * Cbot

    Ex = zeros(Complex128, nLayer)
    Ey = copy(Ex)
    Ez = zeros(Complex128, nLayer-1)
    Hx = copy(Ex)
    Hy = copy(Ex)

    # loop over layers again to compute the fields.
    for j=1:nLayer
        F = Sprod[:,:,j] * MC
        Ex[j] = F[1]
        Ey[j] = F[2]
        Hx[j] = F[3]
        Hy[j] = F[4]

        if j>1
            kp = -iom * xip[j-1]
            km = -iom * xim[j-1]
            QpDerQm = QpNeg[j-1]*QmNegRec[j-1]
            dq = 1-QpDerQm
            dz = 0.5 * zLen[j-1]

            sp = sh(kp*dz)
            sm = sh(km*dz)
            cp = ch(kp*dz)
            cm = ch(km*dz)

            S = zeros(Complex{T}, 4, 4)
            S[1, 1] =  cp-QpDerQm*cm
            S[1, 2] = -QmNegRec[j-1]*(cp-cm)
            S[1, 3] =  QmNegRec[j-1]*(sp/xip[j-1]-sm/xim[j-1])
            S[1, 4] =  sp/xip[j-1]-QpDerQm*sm/xim[j-1]
            S[2, 1] =  QpNeg[j-1]*(cp-cm)
            S[2, 2] = -QpDerQm*cp+cm
            S[2, 3] =  QpDerQm*sp/xip[j-1]-sm/xim[j-1]
            S[2, 4] =  QpNeg[j-1]*(sp/xip[j-1]-sm/xim[j-1])
            S[3, 1] = -QpNeg[j-1]*(xip[j-1]*sp-xim[j-1]*sm)
            S[3, 2] =  QpDerQm*xip[j-1]*sp-xim[j-1]*sm
            S[3, 3] = -QpDerQm*cp+cm
            S[3, 4] = -QpNeg[j-1]*(cp-cm)
            S[4, 1] =  xip[j-1]*sp-QpDerQm*xim[j-1]*sm
            S[4, 2] = -QmNegRec[j-1]*(xip[j-1]*sp-xim[j-1]*sm)
            S[4, 3] =  QmNegRec[j-1]*(cp-cm)
            S[4, 4] =  cp-QpDerQm*cm

            S = S/dq

            # computes E-fields at layer center
            SF = S * F
            Exc = SF[1]
            Eyc = SF[2]
            Ez[j-1] = -(sig[j-1, 5]*Exc + sig[j-1, 6]*Eyc)/sig[j-1, 3]
        end

    end   # nLayer

    if timeFac == "pos"
        conj!(Ex); conj!(Ey); conj!(Ez); conj!(Hx); conj!(Hy)
    end

    if compH
        return Ex, Ey, Ez, Hx, Hy
    end

    return Ex, Ey, Ez

end  # mt1DAniAnalyticField



"""
`mt1DAniAnalyticImp` computes MT responses analytically for 1D general anisotropic model.
The solution is basd on impedance propagation formulas. (see the document by Yuguo Li).
The approach is essentially the same as that of Josef Pek (2002, Computers & Geosciences).

Time dependence: e^{-iwt}.

Input:
    freqs   :: Vector   - frequency array.
    sigma   :: Array    - nLayer*6, conductivity tensor of layered model.
    zNode   :: Vector   - nLayer, depth of top of each layer.

Output:
    Z, Rho, Phs :: Array  - Impedance, apparent resistivity and phase at the top
                            of the first layer, respectively.

"""
function mt1DAniAnalyticImp{T<:Float64}(freqs::Vector{T}, sigma::Array{T}, zNode::Vector{T})

    # check if layer conductivity and depth array have the correct size
    if size(sigma,1) != length(zNode)
        error("The dimension of layer conductivity and depth mismatch!")
    end

    # effective azimuthal anisotropic conductivity
    sigEff = getEffSigma(sigma)

    nFreq  = length(freqs)
    nLayer = size(sigma, 1)

    h = diff(zNode)

    sh = (x::Complex) -> 0.5*(exp(x) - exp(-x))   # sinh()
    ch = (x::Complex) -> 0.5*(exp(x) + exp(-x))   # cosh()

    MU0 = 4*pi*1e-7


    Z   = zeros(Complex128, 4, nFreq)
    Rho = zeros(Float64, 4, nFreq)
    Phs = copy(Rho)

    for iFreq=1:nFreq
        omega = 2 * pi * freqs[iFreq]
        iom = 1im * omega * MU0       # Time dependence: e^{-iwt}

        Zxx=[]; Zxy=[]; Zyx=[]; Zyy=[]
        for j = nLayer:-1:1

            Axx = sigEff[j, 1]
            Ayy = sigEff[j, 2]
            Axy = sigEff[j, 3]

            ada = Axx + Ayy
            add = Axx - Ayy

            dA12 = sqrt( add^2 + 4*Axy^2 )
            if add<0; dA12 = -dA12; end

            A1 = 0.5 * (ada + dA12)
            A2 = 0.5 * (ada - dA12)

            kp = sqrt(-iom) * sqrt(A1);
            km = sqrt(-iom) * sqrt(A2);

            xip = -kp/iom;
            xim = -km/iom;

            QpNeg    = Axy == 0 ? 0 : 2*Axy/(add+dA12)         # -Qp
            QmNegRec = Axy == 0 ? 0 : 0.5*(add-dA12)/Axy       # -1/Qm

            QpDerQm = QpNeg * QmNegRec                         # Q_p/Q_m
            dq = 1-QpDerQm                                     # 1-Q_p/Q_m

            # impedance tensor at the top of the basement
            if j==nLayer
                Zxx =  (1./xip - 1./xim) * QmNegRec/dq
                Zxy =  (1./xip - QpDerQm/xim)/dq
                Zyx = -(1./xim - QpDerQm/xip)/dq
                Zyy =  (1./xip - 1./xim) * QpNeg/dq

                continue
            end

            sp = sh(kp*h[j])
            sm = sh(km*h[j])
            cp = ch(kp*h[j])
            cm = ch(km*h[j])

            S = zeros(Complex{T}, 4, 4)
            S[1, 1] =  cp-QpDerQm*cm
            S[1, 2] = -QmNegRec*(cp-cm)
            S[1, 3] =  QmNegRec*(sp/xip-sm/xim)
            S[1, 4] =  sp/xip-QpDerQm*sm/xim
            S[2, 1] =  QpNeg*(cp-cm)
            S[2, 2] = -QpDerQm*cp+cm
            S[2, 3] =  QpDerQm*sp/xip-sm/xim
            S[2, 4] =  QpNeg*(sp/xip-sm/xim)
            S[3, 1] = -QpNeg*(xip*sp-xim*sm)
            S[3, 2] =  QpDerQm*xip*sp-xim*sm
            S[3, 3] = -QpDerQm*cp+cm
            S[3, 4] = -QpNeg*(cp-cm)
            S[4, 1] =  xip*sp-QpDerQm*xim*sm
            S[4, 2] = -QmNegRec*(xip*sp-xim*sm)
            S[4, 3] =  QmNegRec*(cp-cm)
            S[4, 4] =  cp-QpDerQm*cm

            S = S/dq

            a11 = S[3,1]*Zxx + S[3,2]*Zyx + S[3,3]
            a12 = S[3,1]*Zxy + S[3,2]*Zyy + S[3,4]
            a21 = S[4,1]*Zxx + S[4,2]*Zyx + S[4,3]
            a22 = S[4,1]*Zxy + S[4,2]*Zyy + S[4,4]

            # adds a small value to avoid zero
            deta = a11*a22-a12*a21  + 1.0e-100

            b11 =  a22/deta
            b12 = -a12/deta
            b21 = -a21/deta
            b22 =  a11/deta

            a11 = S[1,1]*Zxx + S[1,2]*Zyx + S[1,3]
            a12 = S[1,1]*Zxy + S[1,2]*Zyy + S[1,4]
            a21 = S[2,1]*Zxx + S[2,2]*Zyx + S[2,3]
            a22 = S[2,1]*Zxy + S[2,2]*Zyy + S[2,4]

            # Z = AB
            Zxx = a11*b11 + a12*b21
            Zxy = a11*b12 + a12*b22
            Zyx = a21*b11 + a22*b21
            Zyy = a21*b12 + a22*b22

        end  # nLayer

        rhoxx = abs(Zxx)^2 / (omega * MU0)
        rhoxy = abs(Zxy)^2 / (omega * MU0)
        rhoyx = abs(Zyx)^2 / (omega * MU0)
        rhoyy = abs(Zyy)^2 / (omega * MU0)

        phsxx = atan2(imag(Zxx), real(Zxx)) * 180/pi
        phsxy = atan2(imag(Zxy), real(Zxy)) * 180/pi
        phsyx = atan2(imag(Zyx), real(Zyx)) * 180/pi
        phsyy = atan2(imag(Zyy), real(Zyy)) * 180/pi

        Z[1:4, iFreq]   = [Zxx Zxy Zyx Zyy]
        Rho[1:4, iFreq] = [rhoxx rhoxy rhoyx rhoyy]
        Phs[1:4, iFreq] = [phsxx phsxy phsyx phsyy]

    end  # nFreq

    return Z, Rho, Phs

end  # mt1DAniFwd



#
# `getEffSigma` computes effective azimuthal anisotropic conductivity
# (see Josef Pek et al., 2002).
function getEffSigma{T<:Real}(sigma::Array{T,2})

    sigEff = zeros(eltype(sigma), size(sigma,1), 3)
    sigEff[:,1] = sigma[:,1] - (sigma[:,5].^2) ./ sigma[:,3]
    sigEff[:,2] = sigma[:,2] - (sigma[:,6].^2) ./ sigma[:,3]
    sigEff[:,3] = sigma[:,4] - (sigma[:,5] .* sigma[:,6]) ./ sigma[:,3]

    return sigEff
end


end # MT1DFwdAni
