"""
`writeModel3DAni` writes 3D general anisotropic EM model parameters.

Input:
    modfile::String       - the name of the model file.
    emMesh::EMTensorMesh  - model parameters and mesh properties.

"""
function writeModel3DAni(modfile::String, emMesh::EMTensorMesh)

    modID = open(modfile,"w")

    @printf(modID,"%-18s %s\n","# Format:", "Model3DAni")
    descrb  = "Model file generated at " * Libc.strftime(time())
    @printf(modID,"%-18s %s\n","# Description:", descrb)

    nx = length(emMesh.xLen)
    ny = length(emMesh.yLen)
    nz = length(emMesh.zLen)

    # cell size along x-direction
    @printf(modID,"%-5s %4d\n","NX:", nx)
    for i = 1:nx
        @printf(modID,"%10.1f", emMesh.xLen[i])
		mod(i,8) == 0 && @printf(modID,"\n")
    end
    mod(nx,8) > 0 && @printf(modID,"\n")

    # cell size along y-direction
    @printf(modID,"%-5s %4d\n","NY:", ny)
    for i = 1:ny
        @printf(modID,"%10.1f", emMesh.yLen[i])
		mod(i,8) == 0 && @printf(modID,"\n")
    end
	mod(ny,8) > 0 && @printf(modID,"\n")

    # air layer
    nAir = 0
    if !isempty(emMesh.airLayer)
        nAir = length(emMesh.airLayer)
        @printf(modID,"%-5s %4d\n","NAIR:", nAir)
        for i = 1:nAir
            @printf(modID,"%10.1f", emMesh.airLayer[i])
			mod(i,8) == 0 && @printf(modID,"\n")
        end
        mod(nAir,8) > 0 && @printf(modID,"\n")
    end


    # cell size along z-direction
    nzEarth = nz-nAir
    @printf(modID,"%-5s %4d\n","NZ:", nzEarth)
    for i = nAir+1:nz
        @printf(modID,"%10.1f", emMesh.zLen[i])
		mod(i-nAir,8) == 0 && @printf(modID,"\n")
    end
	mod(nzEarth,8) != 0 && @printf(modID,"\n")

    # conductivity
    sigma = emMesh.sigma
    if sigma.rotated; rotateTensor3Inv!(sigma); end

    sigxx  = sigma.sig[:,1]
    sigyy  = sigma.sig[:,2]
    sigzz  = sigma.sig[:,3]
    strike = sigma.sig[:,4]
    dip    = sigma.sig[:,5]
    slant  = sigma.sig[:,6]

    if nAir > 0
        airCell = nx * ny * nAir
        sigxx   = sigxx[airCell+1:end]
        sigyy   = sigyy[airCell+1:end]
        sigzz   = sigzz[airCell+1:end]
        strike  = strike[airCell+1:end]
        dip     = dip[airCell+1:end]
        slant   = slant[airCell+1:end]
    end
    sigxx  = reshape(sigxx, nx, ny, nzEarth)
    sigyy  = reshape(sigyy, nx, ny, nzEarth)
    sigzz  = reshape(sigzz, nx, ny, nzEarth)
    strike = reshape(strike, nx, ny, nzEarth)
    dip    = reshape(dip, nx, ny, nzEarth)
    slant  = reshape(slant, nx, ny, nzEarth)
    resType = "Conductivity"
    modType = "Linear"

    @printf(modID,"%-18s %s\n","Resistivity Type:",resType)
    @printf(modID,"%-18s %s\n","Model Type:", modType)

    @printf(modID,"%s\n","Sigma_X:")
    for k = 1:nzEarth
        for j = 1:ny
            for i = 1:nx
                @printf(modID,"%5.2e ",sigxx[i,j,k])
            end
            @printf(modID,"\n")
        end
        k<nzEarth && @printf(modID,"\n")
    end
    @printf(modID,"%s\n","Sigma_Y:")
    for k = 1:nzEarth
        for j = 1:ny
            for i = 1:nx
                @printf(modID,"%5.2e ",sigyy[i,j,k])
            end
            @printf(modID,"\n")
        end
        if k<nzEarth; @printf(modID,"\n"); end
    end
    @printf(modID,"%s\n","Sigma_Z:")
    for k = 1:nzEarth
        for j = 1:ny
            for i = 1:nx
                @printf(modID,"%5.2e ",sigzz[i,j,k])
            end
            @printf(modID,"\n")
        end
        if k<nzEarth; @printf(modID,"\n"); end
    end
    @printf(modID,"%s\n","Sigma_Strike:")
    for k = 1:nzEarth
        for j = 1:ny
            for i = 1:nx
                @printf(modID,"%8.2f ",strike[i,j,k])
            end
            @printf(modID,"\n")
        end
        if k<nzEarth; @printf(modID,"\n"); end
    end
    @printf(modID,"%s\n","Sigma_Dip:")
    for k = 1:nzEarth
        for j = 1:ny
            for i = 1:nx
                @printf(modID,"%8.2f ",dip[i,j,k])
            end
            @printf(modID,"\n")
        end
        if k<nzEarth; @printf(modID,"\n"); end
    end
    @printf(modID,"%s\n","Sigma_Slant:")
    for k = 1:nzEarth
        for j = 1:ny
            for i = 1:nx
                @printf(modID,"%8.2f ",slant[i,j,k])
            end
            @printf(modID,"\n")
        end
        if k<nzEarth; @printf(modID,"\n"); end
    end

    # origin
    origin = copy(emMesh.origin)
    if nAir > 0
        airDep = sum(emMesh.airLayer)
        origin[3] -= airDep
    end
    @printf(modID, "%-15s %4.2e %4.2e %4.2e","Origin (m):",
            origin[1], origin[2], origin[3])

    close(modID)

end
