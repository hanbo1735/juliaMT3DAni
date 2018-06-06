"""
 module `MT3DFwdSolver` defines routines to solve MT 3D anisotropic forward
 modeling problem.

 - by Bo Han, Oct 2016

"""

module MT3DFwdSolver

import YeeGridOperators: CondTensor, EMTensorMesh, setupOperators!
import YeeGridOperators: rotateTensor3!, rotateTensor3Inv!
import MT3DFileIO.DataInfo
import Utilities: ddx, avnc, bilinearInterpMat
import MT1DFwdAni.mt1DAniAnalyticField

import MUMPS: MUMPSfactorization, destroyMUMPS
using LinearSolver

export MT3DFwdComp
# export getBoundaryMT3D, getBoundIndex
# export compRespData, compMT3DResp, compFieldsAtRx


"""
`MT3DFwdComp` is the top function being called to solve the 3D MT forward problem.

Input:
    emMesh    :: EMTensorMesh  - model parameters and mesh properties.
    datInfo   :: DataInfo      - accompanying informations about data.
    linSolver :: String        - which linear solver will be used.
    task      :: String        - determine return `predData` or `fwdResp`.
    saveFac   :: Bool          - whether or not save system matrix decomposition factor (by MUMPS).
    iterTol   :: Real          - error tolerance for iterative solvers. (ignored by direct solvers)

Output:
    predData  :: Vector        - forward data.
    fwdResp   :: Array         - pure forward response.
    eField    :: Array         - solution electric fields (at grid nodes).
    fwdDecomp :: Array         - system matrix decomposition factor (by MUMPS).

"""
function MT3DFwdComp(emMesh::EMTensorMesh, datInfo::DataInfo;
                     linSolver::String="mumps", task::AbstractString="",
                     saveFac::Bool=false, iterTol::Real=1e-6)


    MU0 = 4*pi*1e-7
    muMat = ones(emMesh.nGrid) * MU0

    isempty(emMesh.Curl) && setupOperators!(emMesh)

    # extract things from emMesh
    Curl  = emMesh.Curl
    Grad  = emMesh.Grad
    AveFC = emMesh.AveFC
    AveEC = emMesh.AveEC
    AveNC = emMesh.AveNC
    Vol   = emMesh.Vol
    sigma = emMesh.sigma
    gridSize = emMesh.gridSize
    xLen = emMesh.xLen
    yLen = emMesh.yLen
    zLen = emMesh.zLen

    !sigma.rotated && rotateTensor3!(sigma)

    sigxx = sigma.sig[:,1]
    sigyy = sigma.sig[:,2]
    sigzz = sigma.sig[:,3]
    sigxy = sigma.sig[:,4]
    sigxz = sigma.sig[:,5]
    sigyz = sigma.sig[:,6]

    nx, ny, nz = gridSize
    nNode = (nx+1)*(ny+1)*(nz+1)
    nEx = nx*(ny+1)*(nz+1)
    nEy = (nx+1)*ny*(nz+1)
    nEz = (nx+1)*(ny+1)*nz


    # setup mass matrices
    MmuF = AveFC' * (Vol * (1 ./ muMat))
    MmuF = spdiagm(MmuF)

    MmuN = AveNC' * (Vol * (1 ./ muMat))
    MmuN = spdiagm(MmuN)

    # coefficients arising from diagonal anisotropy
    AveX = AveEC[:,1:nEx]
    AveY = AveEC[:,nEx+1:nEx+nEy]
    AveZ = AveEC[:,nEx+nEy+1:end]
    Msig = vcat(AveX'*(Vol*(sigxx)), AveY'*(Vol*(sigyy)), AveZ'*(Vol*(sigzz)))
    Msig = spdiagm(Msig)

    # coefficients arising from off-diagonal anisotropy
    MsigAni = getAnisoOperator(gridSize, Vol, sigxy, sigxz, sigyz)


    # frequencies
    freqs = datInfo.freqs
    omega = 2 * pi * freqs
    nFreq = length(omega)

    fwdDecomp = []

    if uppercase(linSolver) == "MUMPS"
        lsFlag = 1

        # save MUMPS matrix decomposition facor
        saveFac && ( fwdDecomp = Array{MUMPSfactorization}(nFreq) )

    elseif uppercase(linSolver) == "PARDISO"
        lsFlag = 2

    else
        # iteration log file
        lsFlag = 3
        iterLogID = open("iter.log", "w")

    end


    # solved electric and magnetic fields of two polarization modes
    (nF, nE) = size(Curl)
    eField = Array{Complex{Float64}}(nE, 2, nFreq)
    bField = Array{Complex{Float64}}(nF, 2, nFreq)

    # edge divergence operator
    # DivE  = -Grad'

    # Laplacian operator:   ∇×μ^{-1}∇× - ∇μ^{-1}∇
    # Lap = (Curl' * MmuF) * Curl - Grad * MmuN * DivE

    dCurl = (Curl' * MmuF) * Curl

    ii, io = getBoundIndex(gridSize::Vector{Int})
    rAii = dCurl[ii, ii]
    rAio = dCurl[ii, io]
    iAii = (Msig + MsigAni)[ii, ii]
    iAio = (Msig + MsigAni)[ii, io]

    # get the inner node gradient. used by Aphi equation
    idx3D  = reshape(collect(1:nNode), nx+1, ny+1, nz+1)
    iiND   = idx3D[2:end-1, 2:end-1, 2:end-1][:]
    Gradii = Grad[ii, iiND]

    STBa = Gradii * MmuN[iiND, iiND] * Gradii'

    # compute boundary fields for two polarization modes
    bcXY = getBoundaryMT3D(freqs, xLen, yLen, zLen, sigma.sig, 1)
    bcYX = getBoundaryMT3D(freqs, xLen, yLen, zLen, sigma.sig, 2)

    eField[io, 1, :] = bcXY[:, :]
    eField[io, 2, :] = bcYX[:, :]


    # loop over frequency
    for i = 1:nFreq
        iome = 1im * omega[i]

        # coefficients for inner electric fields
        Aii = rAii + iome * iAii

        # coefficients for boundary electric fields
        Aio = rAio + iome * iAio

        # rhs contains two polarization modes
        rhs = -Aio * hcat(bcXY[:,i], bcYX[:,i])

        # solve the linear system
        println("Solving the forward problem of the $i / $nFreq th frequency with $(uppercase(linSolver)) ...")
        if lsFlag == 1      # MUMPS
            if saveFac
                @time etmp, fwdDecomp[i] = mumpsSolver(Aii, rhs, saveFac=true)
            else
                @time etmp = mumpsSolver(Aii, rhs)
            end

        elseif lsFlag == 2  # PARDISO
            # not developed yet
            @time etmp = pardiSolver(Aii, rhs)

        elseif lsFlag == 3  # iterative solver
            # set up preconditioner

            # Jacobi and SSOR preconditioner
            # D  = diag(Aii)
            # PC(x) = D.\x
            # PC(x) = tril(Aii)\(D.*(triu(Aii)\x))


            # SSOR preconditioner based on A-phi system
            Aap = [Aii + STBa               iome * iAii * Gradii;
                   iome * Gradii' * iAii    iome * Gradii' * iAii * Gradii]

            Map(x) = tril(Aap)\(diag(Aap).*(triu(Aap)\x))
            P1     = [speye(size(Aii,2));  Gradii']
            P2     = [speye(size(Aii,1))   Gradii]
            PC(x)  = P2*(Map(P1*x))


            @printf(iterLogID,"%s %3d\n","Freq No.:", i)
            @time etmp = iterSolver(Aii, rhs, PC, linSolver, iterLogID, iterTol)

        end

        eField[ii,:,i] = etmp

        # lsFlag == 1 && destroyMUMPS(Ainv)

        tmp = view(eField, 1:nE, 1:2, i)
        bField[:, :, i] = -1 / (iome) * (Curl * tmp)

    end

    lsFlag == 3 && close(iterLogID)

    # Pure forward problem
    if uppercase(task) == "PF"
        fwdResp = compRespData(emMesh, datInfo, eField, bField, task)
        return fwdResp, eField
    end

    # compute predicted data at receivers
    predData = compRespData(emMesh, datInfo, eField, bField)

    if lsFlag == 1 && saveFac   # MUMPS
        return predData, eField, fwdDecomp
    else
        return predData, eField
    end

end


include("boundaryFields.jl")
include("dataFunc.jl")
include("getAnisoOperator.jl")

end # MT3DFwdSolver
