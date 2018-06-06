"""
`readModel3DAni` reads 3D general anisotropic EM model parameters.

Input:
    modfile::String       - the name of the model file.

Output:
    emMesh::EMTensorMesh  - model parameters and mesh properties.

"""
function readModel3DAni(modfile::String)

    if isfile(modfile)
        fid = open(modfile,"r")
    else
        error("$(modfile) does not exist, please try again!")
    end

    # declare variables to be output
    xLen = []
    yLen = []
    zLen = []
    nx   = []
    ny   = []
    nz   = []
    nAir = 0
    airLayer = zeros(0)
    sigX = []
    sigY = []
    sigZ = []
    strike = []
    dip    = []
    slant  = []
    origin = []

    while !eof(fid)

        cline = strip(readline(fid))

        # ignore empty lines and comments (preceded with #)
        isempty(cline)  && continue
        cline[1] == '#' && continue

        # cell size along x-direction
        if contains(cline, "NX")
            tmp = split(cline)
            nx  = parse(Int, tmp[end])
            nd  = 0
            xLen = zeros(nx)

            while nd < nx
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    xLen[nd] = float(cline[i])
                end
            end

        # cell size along y-direction
        elseif contains(cline, "NY")
            tmp = split(cline)
            ny  = parse(Int, tmp[end])
            nd  = 0
            yLen = zeros(ny)

            while nd < ny
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    yLen[nd] = float(cline[i])
                end
            end

        # cell size along z-direction
        elseif contains(cline, "NZ")
            tmp = split(cline)
            nz  = parse(Int, tmp[end])
            nd  = 0
            zLen = zeros(nz)

            while nd < nz
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    zLen[nd] = float(cline[i])
                end
            end

        # air layer size along z-direction
        elseif contains(cline, "NAIR")
            tmp  = split(cline)
            nAir = parse(Int, tmp[end])
            nd   = 0
            airLayer = zeros(nAir)

            while nd < nAir
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    airLayer[nd] = float(cline[i])
                end
            end

        # resistivity type: resistivity or conductivity
        elseif contains(cline, "Resistivity Type")
            tmp = split(cline)
            resType = lowercase(tmp[end])

        # model type: linear or logorithmic
        elseif contains(cline, "Model Type")
            tmp = split(cline)
            modType = lowercase(tmp[end])
            nBlock = nx * ny * nz

            sigX = zeros(nBlock)
            nd = 0
            cline = readline(fid)
            while nd < nBlock
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    sigX[nd] = float(cline[i])
                end
            end

            sigY = zeros(nBlock)
            nd = 0
            cline = readline(fid)
            while nd < nBlock
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    sigY[nd] = float(cline[i])
                end
            end

            sigZ = zeros(nBlock)
            nd = 0
            cline = readline(fid)
            while nd < nBlock
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    sigZ[nd] = float(cline[i])
                end
            end

            strike = zeros(nBlock)
            nd = 0
            cline = readline(fid)
            while nd < nBlock
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    strike[nd] = float(cline[i])
                end
            end

            dip = zeros(nBlock)
            nd = 0
            cline = readline(fid)
            while nd < nBlock
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    dip[nd] = float(cline[i])
                end
            end

            slant = zeros(nBlock)
            nd = 0
            cline = readline(fid)
            while nd < nBlock
                cline = strip(readline(fid))
                cline = split(cline)
                num = length(cline)
                for i = 1:num
                    nd = nd + 1
                    slant[nd] = float(cline[i])
                end
            end

            # convert to linear conductivity
            if modType == "log"
                sigX = exp(sigX)
                sigY = exp(sigY)
                sigZ = exp(sigZ)
            end
            if resType == "resistivity"
                sigX = 1 ./ sigX
                sigY = 1 ./ sigY
                sigZ = 1 ./ sigZ
            end

        # origin
        elseif contains(cline, "Origin")
            tmp = split(cline)
            origin = zeros(3)
            origin[1] = float(tmp[end-2])
            origin[2] = float(tmp[end-1])
            origin[3] = float(tmp[end])

        end

    end

    close(fid)

    # append air layers, note that airLayer is from down to up
    airThickness = sum(airLayer)
    zLen = vcat(flipdim(airLayer, 1), zLen)
    origin[3] = origin[3] + airThickness
    sigAir = 1e-8
    airVec = ones(nx*ny*nAir) * sigAir
    sigX   = vcat(airVec, sigX)
    sigY   = vcat(airVec, sigY)
    sigZ   = vcat(airVec, sigZ)

    ztmp   = zeros(nx*ny*nAir)
    strike = vcat(ztmp, strike)
    dip    = vcat(ztmp, dip)
    slant  = vcat(ztmp, slant)

    nz = length(zLen)

    sigFull = hcat(sigX, sigY, sigZ, strike, dip, slant)
    sigma = CondTensor(sigFull, false)

    gridSize = [nx, ny, nz]
    nGrid = nx * ny * nz
    empMat = spzeros(0,0)

    emMesh = EMTensorMesh(xLen, yLen, zLen, airLayer, gridSize,
                          origin, nGrid, sigma, empMat, empMat,
                          empMat, empMat, empMat, empMat, empMat)

    return emMesh

end
