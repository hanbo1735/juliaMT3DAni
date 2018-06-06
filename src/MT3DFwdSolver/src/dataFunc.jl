"""
`compRespData` computes predicted forward data at receiver locations.

Input:
    emMesh  :: EMTensorMesh  - model parameters and mesh properties.
    datInfo :: DataInfo      - accompanying informations about data.
    eField  :: Array         - solution electric fields (at grid nodes).
    bField  :: Array         - solution magnetic fields (at grid nodes).
    task    :: String        - determine return `predData` or `fwdResp`.

Output:
    predData :: Vector    - forward data.
    fwdResp  :: Array     - pure forward response.

"""
function compRespData{T<:Complex}(emMesh::EMTensorMesh, datInfo::DataInfo,
                                  eField::Array{T,3}, bField::Array{T,3},
                                  task::AbstractString="")

    xLen = emMesh.xLen
    yLen = emMesh.yLen
    zLen = emMesh.zLen

    rxID     = datInfo.rxID
    dcID     = datInfo.dcID
    freqID   = datInfo.freqID
    dataType = datInfo.dataType
    dataComp = datInfo.dataComp

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    rxLoc = datInfo.rxLoc[:, 1:2]
    xNode = cumsum([0; xLen]) - emMesh.origin[1]
    yNode = cumsum([0; yLen]) - emMesh.origin[2]
    nAir  = length(emMesh.airLayer)
    zLen1 = emMesh.zLen[nAir+1]
    sigma = emMesh.sigma.sig

    sigma1 = zeros(Float64, nx, ny, 6)
    sigma1[:, :, 1] = reshape(sigma[:,1], nx, ny, nz)[:, :, nAir+1]
    sigma1[:, :, 2] = reshape(sigma[:,2], nx, ny, nz)[:, :, nAir+1]
    sigma1[:, :, 3] = reshape(sigma[:,3], nx, ny, nz)[:, :, nAir+1]
    sigma1[:, :, 4] = reshape(sigma[:,4], nx, ny, nz)[:, :, nAir+1]
    sigma1[:, :, 5] = reshape(sigma[:,5], nx, ny, nz)[:, :, nAir+1]
    sigma1[:, :, 6] = reshape(sigma[:,6], nx, ny, nz)[:, :, nAir+1]

    nEx = nx * (ny+1) * (nz+1)
    nEy = (nx+1) * ny * (nz+1)
    nEz = (nx+1) * (ny+1) * nz

    nBx = (nx+1) * ny * nz
    nBy = nx * (ny+1) * nz
    nBz = nx * ny * (nz+1)

    nRx   = size(rxLoc, 1)
    nFreq = length(datInfo.freqs)
    if contains(dataType, "Impedance")
        if contains(dataType, "Tipper")
            fwdResp = zeros(Complex{Float64}, nRx, 6, nFreq)
        else
            fwdResp = zeros(Complex{Float64}, nRx, 4, nFreq)
        end
    elseif contains(dataType, "Rho_Phs")
        if contains(dataType, "Tipper")
            fwdResp = zeros(Float64, nRx, 12, nFreq)
        else
            fwdResp = zeros(Float64, nRx, 8, nFreq)
        end
    end

    for i=1:nFreq

        # 1st interpolate EM fields at grid nodes to receiver locations.
        nPol = 2
        EHr = zeros(Complex{Float64}, nRx, 5, nPol)

        for j=1:nPol

            ex = reshape(eField[1:nEx, j, i], nx, ny+1, nz+1)
            ey = reshape(eField[nEx+1:nEx+nEy, j, i], nx+1, ny, nz+1)
            ez = reshape(eField[nEx+nEy+1:nEx+nEy+nEz, j, i], nx+1, ny+1, nz)
            bx = reshape(bField[1:nBx, j, i], nx+1, ny, nz)
            by = reshape(bField[nBx+1:nBx+nBy, j, i], nx, ny+1, nz)
            bz = reshape(bField[nBx+nBy+1:end, j, i], nx, ny, nz+1)

            ex12 = ex[:, :, nAir+1:nAir+2]
            ey12 = ey[:, :, nAir+1:nAir+2]
            bz12 = bz[:, :, nAir+1:nAir+2]
            ez1  = ez[:, :, nAir+1]
            bx1  = bx[:, :, nAir+1]
            by1  = by[:, :, nAir+1]

            EHr[:, :, j] = compFieldsAtRx(rxLoc, xNode, yNode, zLen1, sigma1,
                                          ex12, ey12, ez1, bx1, by1, bz12)
        end # nPol

        # 2nd compute MT responses from EM fields.
        freq = datInfo.freqs[i]
        fwdResp[:, :, i] = compMT3DResp(freq, EHr, dataType)

    end  # nFreq

    # if it's a pure forward problem, then we finished here.
    uppercase(task) == "PF" && return fwdResp


    # it's not a pure forward problem by default, then we need to select
    # the predicted data
    nData = length(rxID)
    if contains(dataType, "Impedance")
        predData = zeros(Complex{Float64}, nData)
    elseif contains(dataType, "Rho_Phs")
        predData = zeros(Float64, nData)
    end

    p = 1
    for i = 1:nFreq

        indF  = find(freqID .== i)

        # check if current frequency exists
        isempty(indF) && continue

        subRxID = rxID[indF]
        subDcID = dcID[indF]
        nd  = length(subRxID)

        for j = 1:nd

            iRx = subRxID[j]
            iDt = subDcID[j]
            dt  = dataComp[iDt]

            if contains(dataType, "Impedance")

                if dt == "ZXX"
                    predData[p:p] = fwdResp[iRx, 1, i]
                elseif dt == "ZXY"
                    predData[p:p] = fwdResp[iRx, 2, i]
                elseif dt == "ZYX"
                    predData[p:p] = fwdResp[iRx, 3, i]
                elseif dt == "ZYY"
                    predData[p:p] = fwdResp[iRx, 4, i]
                end

                if contains(dataType, "Tipper")
                    if dt == "TZX"
                        predData[p:p] = fwdResp[iRx, 5, i]
                    elseif dt == "TZY"
                        predData[p:p] = fwdResp[iRx, 6, i]
                    end
                end

            elseif contains(dataType, "Rho_Phs")

                if dt == "RhoXX"
                    predData[p:p] = fwdResp[iRx, 1, i]
                elseif dt == "PhsXX"
                    predData[p:p] = fwdResp[iRx, 2, i]
                elseif dt == "RhoXY"
                    predData[p:p] = fwdResp[iRx, 3, i]
                elseif dt == "PhsXY"
                    predData[p:p] = fwdResp[iRx, 4, i]
                elseif dt == "RhoYX"
                    predData[p:p] = fwdResp[iRx, 5, i]
                elseif dt == "PhsYX"
                    predData[p:p] = fwdResp[iRx, 6, i]
                elseif dt == "RhoYY"
                    predData[p:p] = fwdResp[iRx, 7, i]
                elseif dt == "PhsYY"
                    predData[p:p] = fwdResp[iRx, 8, i]
                end

                if contains(dataType, "Tipper")
                    if dt == "RealTZX"
                        predData[p:p] = fwdResp[iRx, 9, i]
                    elseif dt == "ImagTZX"
                        predData[p:p] = fwdResp[iRx, 10, i]
                    elseif dt == "RealTZY"
                        predData[p:p] = fwdResp[iRx, 11, i]
                    elseif dt == "ImagTZY"
                        predData[p:p] = fwdResp[iRx, 12, i]
                    end
                end

            else
                error("Data type $(dataType) is not supported!")

            end

            p += 1

        end # j = nd

    end # i = nFreq

    return predData

end



"""
`compFieldsAtRx` computes EM fields at receiver locations from fields at grid nodes
for a single frequency with a single polarization mode.

Input:
    rxLoc  :: Array     - receiver locations.
    xNode  :: Vector    - node x-coordinates.
    yNode  :: Vector    - node y-coordinates.
    zLen1  :: Float64   - z-cell size of the receiver layer.
    sigma1 :: Array     - cell conductivity of the receiver layer.
    Ex, Ey, Bz :: Array{Complex128, 3}    - grid fields of the receiver layer.
    Ez, Bx, By :: Array{Complex128, 2}    - grid fields of the receiver layer.

Output:
    EHr :: Array{Complex128,nRx,5}   - Ex, Ey, Hx, Hy, Hz at the receiver locations.

"""
function compFieldsAtRx{T1<:Float64, T2<:Complex128}(rxLoc::Array{T1,2}, xNode::Vector{T1},
                                         yNode::Vector{T1}, zLen1::T1, sigma1::Array{T1,3},
                                         Ex::Array{T2,3}, Ey::Array{T2,3}, Ez::Array{T2,2},
                                         Bx::Array{T2,2}, By::Array{T2,2}, Bz::Array{T2,3})


    MU0 = 4 * pi * 1e-7

    xLen = diff(xNode)
    yLen = diff(yNode)

    # cell center coordiantes
    xCen = xNode[1:end-1] + xLen / 2
    yCen = yNode[1:end-1] + yLen / 2

    nx = length(xLen)
    ny = length(yLen)

    sig1xx = sigma1[:,:,1]
    sig1yy = sigma1[:,:,2]
    sig1zz = sigma1[:,:,3]
    sig1xy = sigma1[:,:,4]
    sig1xz = sigma1[:,:,5]
    sig1yz = sigma1[:,:,6]

    # First compute Hy at the receiver layer (earth surface or seafloor).
    Hy0 = zeros(Complex128, nx, ny+1)

    for ix=1:nx
        Bz0 = vec(Bz[ix, :, 1])
        Bz1 = vec(Bz[ix, :, 2])

        Ex0 = vec(Ex[ix, 2:end-1, 1])
        Ex1 = vec(Ex[ix, 2:end-1, 2])

        Ey0 = 0.5 * ( vec(Ey[ix, :, 1]) + vec(Ey[ix+1, :, 1]) )
        Ey1 = 0.5 * ( vec(Ey[ix, :, 2]) + vec(Ey[ix+1, :, 2]) )

        EzH = 0.5 * ( vec(Ez[ix, 2:end-1]) + vec(Ez[ix+1, 2:end-1]) )

        # quarter Hz (1/4), with length ny
        HzQ = (0.75 * Bz0 + 0.25 * Bz1) / MU0

        # half Hy (1/2), with length ny-1
        # More strictly, an average mu should be used here.
        HyH = vec(By[ix, 2:end-1]) / MU0

        # quarter Ex (1/4), with length ny-1
        ExQ = 0.75 * Ex0 + 0.25 * Ex1

        # quarter Ey collocated with ExQ
        Ey0t = (avnc(ny-1)*(Ey0./yLen)) .* (yLen[1:end-1].*yLen[2:end]) ./ (avnc(ny-1)*yLen)
        Ey1t = (avnc(ny-1)*(Ey1./yLen)) .* (yLen[1:end-1].*yLen[2:end]) ./ (avnc(ny-1)*yLen)
        EyQ  = 0.75 * Ey0t + 0.25 * Ey1t

        # quarter Ez collocated with ExQ, assume Ez=0 at the earth surface
        EzQ = 0.5 * EzH


        # average conductiviy at vertical edge
        sigtxx = vec(sig1xx[ix,:])
        sigtxy = vec(sig1xy[ix,:])
        sigtxz = vec(sig1xz[ix,:])

        sig1xxV = (avnc(ny-1)*(sigtxx.*yLen)) ./ (avnc(ny-1)*yLen)
        sig1xyV = (avnc(ny-1)*(sigtxy.*yLen)) ./ (avnc(ny-1)*yLen)
        sig1xzV = (avnc(ny-1)*(sigtxz.*yLen)) ./ (avnc(ny-1)*yLen)

        JxQ = sig1xxV.*ExQ + sig1xyV.*EyQ + sig1xzV.*EzQ

        # dHz/dy
        dHzQ = (ddx(ny-1)*HzQ) ./ (avnc(ny-1)*yLen)

        # Ampre's theorem: dHz/dy - dHy/dz = JxQ
        # where dHy/dz = (HyH-Hy0)/(0.5*zLen1).
        Hy0[ix, 2:end-1] = HyH - (dHzQ - JxQ)*(0.5*zLen1)
        Hy0[ix, 1]   = Hy0[ix, 2]
        Hy0[ix, end] = Hy0[ix, end-1]
    end  # nx


    # Second compute Hx at the receiver layer (earth surface or seafloor).
    Hx0 = zeros(Complex128, nx+1, ny)

    for iy=1:ny
        Bz0 = vec(Bz[:, iy, 1])
        Bz1 = vec(Bz[:, iy, 2])

        Ey0 = vec(Ey[2:end-1, iy, 1])
        Ey1 = vec(Ey[2:end-1, iy, 2])

        Ex0 = 0.5 * ( vec(Ex[:, iy, 1]) + vec(Ex[:, iy+1, 1]) )
        Ex1 = 0.5 * ( vec(Ex[:, iy, 2]) + vec(Ex[:, iy+1, 2]) )

        EzH = 0.5 * ( vec(Ez[2:end-1, iy]) + vec(Ez[2:end-1, iy+1]) )

        # quarter Hz (1/4), with length nx
        HzQ = (0.75 * Bz0 + 0.25 * Bz1) / MU0

        # half Hx (1/2), with length nx-1
        # More strictly, an average mu should be used here.
        HxH = vec(Bx[2:end-1, iy]) / MU0

        # quarter Ey (1/4), with length nx-1
        EyQ = 0.75 * Ey0 + 0.25 * Ey1

        # quarter Ex collocated with EyQ
        Ex0t = (avnc(nx-1)*(Ex0./xLen)) .* (xLen[1:end-1].*xLen[2:end]) ./ (avnc(nx-1)*xLen)
        Ex1t = (avnc(nx-1)*(Ex1./xLen)) .* (xLen[1:end-1].*xLen[2:end]) ./ (avnc(nx-1)*xLen)
        ExQ  = 0.75 * Ex0t + 0.25 * Ex1t

        # quarter Ez collocated with EyQ, assume Ez=0 at the earth surface
        EzQ = 0.5 * EzH

        # average conductiviy at vertical edge
        sigtyy = vec(sig1yy[:,iy])
        sigtxy = vec(sig1xy[:,iy])
        sigtyz = vec(sig1yz[:,iy])

        sig1yyV = (avnc(nx-1)*(sigtyy.*xLen)) ./ (avnc(nx-1)*xLen)
        sig1xyV = (avnc(nx-1)*(sigtxy.*xLen)) ./ (avnc(nx-1)*xLen)
        sig1yzV = (avnc(nx-1)*(sigtyz.*xLen)) ./ (avnc(nx-1)*xLen)

        JyQ = sig1xyV.*ExQ + sig1yyV.*EyQ + sig1yzV.*EzQ

        # dHz/dx
        dHzQ = (ddx(nx-1)*HzQ) ./ (avnc(nx-1)*xLen)

        # Ampre's theorem: dHx/dz - dHz/dx = JyQ,
        # where dHx/dz = (HxH-Hx0)/(0.5*zLen1).
        Hx0[2:end-1, iy] = HxH - (dHzQ + JyQ)*(0.5*zLen1)
        Hx0[1, iy]   = Hx0[2, iy]
        Hx0[end, iy] = Hx0[end-1, iy]
    end  # ny


    # Third interpolate fields to receiver locations (using bilinear interpolation)
    # pre-defined fields at receiver locations
    nRx = size(rxLoc,1)
    Exr = zeros(Complex128, nRx)
    Eyr = copy(Exr)
    Hxr = copy(Exr)
    Hyr = copy(Exr)
    Hzr = copy(Exr)


    itpMat1 = bilinearInterpMat(rxLoc, xCen, yNode)    # for Ex & Hy
    itpMat2 = bilinearInterpMat(rxLoc, xNode, yCen)    # for Ey & Hx
    itpMat3 = bilinearInterpMat(rxLoc, xCen, yCen)     # for Hz

    Exr = itpMat1' * vec(Ex[:, :, 1])
    Eyr = itpMat2' * vec(Ey[:, :, 1])
    Hxr = itpMat2' * vec(Hx0)
    Hyr = itpMat1' * vec(Hy0)
    Hzr = itpMat3' * vec(Bz[:, :, 1]) / MU0

    return hcat(Exr, Eyr, Hxr, Hyr, Hzr)
end



"""
`compMT3DResp` computes 3D MT responses from EM fields.

Input:
    freq :: Float64      - a single frequency value.
    EHr  :: Array        - EM fields at the receiver locations.
    dataType :: String   - the type of responses (impedance or rho & phase).

Output:
    resp :: Array        - forward responses.
"""
function compMT3DResp{T<:Complex128}(freq::Float64, EHr::Array{T, 3}, dataType::AbstractString)

    MU0 = 4 * pi * 1e-7
    omega = 2 * pi * freq

    ex1 = EHr[:, 1, 1]
    ey1 = EHr[:, 2, 1]
    hx1 = EHr[:, 3, 1]
    hy1 = EHr[:, 4, 1]
    hz1 = EHr[:, 5, 1]

    ex2 = EHr[:, 1, 2]
    ey2 = EHr[:, 2, 2]
    hx2 = EHr[:, 3, 2]
    hy2 = EHr[:, 4, 2]
    hz2 = EHr[:, 5, 2]


    hinv = hx1 .* hy2 - hy1 .* hx2

    Zxx = (ex1 .* hy2 - ex2 .* hy1) ./ hinv
    Zxy = (ex2 .* hx1 - ex1 .* hx2) ./ hinv
    Zyx = (ey1 .* hy2 - ey2 .* hy1) ./ hinv
    Zyy = (ey2 .* hx1 - ey1 .* hx2) ./ hinv

    # Tipper
    #   [Hz1 Hz2] = [Tzx Tzy] * [Hx1 Hx2]
    #                           [Hy1 Hy2]

    Tzx = (hy2 .* hz1 - hy1 .* hz2) ./ hinv
    Tzy = (hx1 .* hz2 - hx2 .* hz1) ./ hinv

    if dataType == "Impedance"
        resp = hcat(Zxx, Zxy, Zyx, Zyy)

    elseif dataType == "Impedance_Tipper"
        resp = hcat(Zxx, Zxy, Zyx, Zyy, Tzx, Tzy)

    elseif contains(dataType, "Rho_Phs")
        rhoxx = (abs.(Zxx)).^2 / (omega*MU0)
        rhoxy = (abs.(Zxy)).^2 / (omega*MU0)
        rhoyx = (abs.(Zyx)).^2 / (omega*MU0)
        rhoyy = (abs.(Zyy)).^2 / (omega*MU0)

        phsxx = atan2.(imag(Zxx), real(Zxx)) * 180/pi
        phsxy = atan2.(imag(Zxy), real(Zxy)) * 180/pi
        phsyx = atan2.(imag(Zyx), real(Zyx)) * 180/pi
        phsyy = atan2.(imag(Zyy), real(Zyy)) * 180/pi

        if contains(dataType, "Tipper")
            resp = hcat(rhoxx, phsxx, rhoxy, phsxy, rhoyx, phsyx, rhoyy, phsyy,
                        real(Tzx), imag(Tzx), real(Tzy), imag(Tzy))
        else
            resp = hcat(rhoxx, phsxx, rhoxy, phsxy, rhoyx, phsyx, rhoyy, phsyy)
        end

    end

    return resp
end
