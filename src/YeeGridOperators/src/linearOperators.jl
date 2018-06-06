
# getFaceDivergence constructs face divergence mapping from
# face center to cell center
function getFaceDivergence{T<:Real}(xLen::Vector{T}, yLen::Vector{T},
                                    zLen::Vector{T})
    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    Dx = kron3(speye(nz), speye(ny), ddx(nx))
    Dy = kron3(speye(nz), ddx(ny), speye(nx))
    Dz = kron3(ddx(nz), speye(ny), speye(nx))

    faceDivMat = [Dx Dy Dz]

    faceArea = meshGeoFace(xLen, yLen, zLen)
    volMat   = meshGeoVolume(xLen, yLen, zLen)
    volMat   = sdInv(volMat)

    faceDivMat  = volMat * faceDivMat * faceArea

    return faceDivMat

end

# getNodalGradient constructs nodal gradient mapping from node to
# edge center
function getNodalGradient{T<:Real}(xLen::Vector{T}, yLen::Vector{T},
                                   zLen::Vector{T})

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    Gx = kron3(speye(nz+1), speye(ny+1), ddx(nx))
    Gy = kron3(speye(nz+1), ddx(ny), speye(nx+1))
    Gz = kron3(ddx(nz), speye(ny+1), speye(nx+1))

    GradMat = [Gx; Gy; Gz]

    edgeMat = meshGeoEdge(xLen, yLen, zLen)
    edgeMat = sdInv(edgeMat)

    GradMat = edgeMat * GradMat

    return GradMat

end

# getCellGradient constructs cell gradient mapping from cell-center to
# cell face
function getCellGradient{T<:Real}(xLen::Vector{T}, yLen::Vector{T},
                                  zLen::Vector{T})

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    Gx = kron3(speye(nz), speye(ny), ddxC2N(nx))
    Gy = kron3(speye(nz), ddxC2N(ny), speye(nx))
    Gz = kron3(ddxC2N(nz), speye(ny), speye(nx))

    cellGradMat = [Gx; Gy; Gz]

    # get face areas and cell volumes
    faceMat = meshGeoFace(xLen, yLen, zLen)
    volMat  = meshGeoVolume(xLen, yLen, zLen)
    AvCF    = aveCell2Face([nx, ny, nz])
    volMat  = AvCF * diag(volMat);

    disMat   = diag(faceMat) ./ volMat
    cellGradMat = sdiag(disMat) * cellGradMat

    nxFace = size(Gx,1)
    nyFace = size(Gy,1)

    xcGrad = cellGradMat[1:nxFace, :]
    ycGrad = cellGradMat[nxFace+1:nxFace + nyFace, :]
    zcGrad = cellGradMat[nxFace + nyFace + 1 : end, :]

    return cellGradMat, xcGrad, ycGrad, zcGrad

end

# getCellGradientBC constructs cell gradient with boundary conditions mapping
# from cell-center to cell face
function getCellGradientBC{T<:Real}(xLen::Vector{T}, yLen::Vector{T},
                                    zLen::Vector{T},BC::Vector{String})

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    #BC = "neumann"

    Dx = kron3(speye(nz), speye(ny), ddxCellGradBC(nx,BC[1]))
    Dy = kron3(speye(nz), ddxCellGradBC(ny,BC[2]), speye(nx))
    Dz = kron3(ddxCellGradBC(nz,BC[3]), speye(ny), speye(nx))

    cellGradMat = [Dx; Dy; Dz]

    # get face areas and cell volumes
    faceMat = meshGeoFace(xLen, yLen, zLen)
    volMat  = meshGeoVolume(xLen, yLen, zLen)
    AvCF    = aveCell2Face([nx, ny, nz])
    volMat  = AvCF * diag(volMat)

    disMat   = diag(faceMat) ./ volMat
    cellGradMat = sdiag(disMat) * cellGradMat

    return cellGradMat

end

# getEdgeCurl constructs edge curl mapping from cell-edge to cell-face
function getEdgeCurl{T<:Real}(xLen::Vector{T}, yLen::Vector{T},
                              zLen::Vector{T})

    nx = length(xLen)
    ny = length(yLen)
    nz = length(zLen)

    nxFace = (nx+1) * ny * nz
    nyFace = nx * (ny+1) * nz
    nzFace = nx * ny * (nz+1)
    nxEdge = nx * (ny+1) * (nz+1)
    nyEdge = (nx+1) * ny * (nz+1)
    nzEdge = (nx+1) * (ny+1) * nz

    ## The curl operator from edges to faces
    Dyz = kron3(-ddx(nz), speye(ny), speye(nx+1))
    Dzy = kron3(speye(nz), -ddx(ny), speye(nx+1))

    Dxz = kron3(-ddx(nz), speye(ny+1), speye(nx))
    Dzx = kron3(speye(nz), speye(ny+1), -ddx(nx))

    Dxy = kron3(speye(nz+1), -ddx(ny), speye(nx))
    Dyx = kron3(speye(nz+1), speye(ny), -ddx(nx))

    ExCurl = [spzeros(nxFace, nxEdge) Dyz  -Dzy]
    EyCurl = [-Dxz spzeros(nyFace, nyEdge)  Dzx]
    EzCurl = [Dxy -Dyx  spzeros(nzFace, nzEdge)]

    EdgeCurlMat = [ExCurl; EyCurl; EzCurl]

    # face area and edge size
    faceMat = meshGeoFace(xLen, yLen, zLen)
    edgeMat = meshGeoEdge(xLen, yLen, zLen)

    faceMat = sdInv(faceMat)

    EdgeCurlMat = faceMat * EdgeCurlMat * edgeMat

    return EdgeCurlMat

end


#----------------------------------------------------------
# Averaging operators over mesh
#----------------------------------------------------------

# aveFace2Cell gets averaging mapping from face-center to cell-center
function aveFace2Cell(n::Vector{Int})

    AvFCx = kron3(speye(n[3]), speye(n[2]), av(n[1]))
    AvFCy = kron3(speye(n[3]), av(n[2]), speye(n[1]))
    AvFCz = kron3(av(n[3]), speye(n[2]), speye(n[1]))

    AveFCMat = [AvFCx AvFCy AvFCz]

    return AveFCMat

end

# aveCell2Face gets averaging mapping from cell-center to face-center
function aveCell2Face(n::Vector{Int})

    AvCFx = kron3(speye(n[3]), speye(n[2]), avcn(n[1]))
    AvCFy = kron3(speye(n[3]), avcn(n[2]), speye(n[1]))
    AvCFz = kron3(avcn(n[3]), speye(n[2]), speye(n[1]))

    AveCFMat = [AvCFx; AvCFy; AvCFz]

    return AveCFMat

end

# aveCell2Edge gets averaging mapping from cell-center to edge
function aveCell2Edge(n::Vector{Int})

    Avcex = avcn(n[1])
    Avcey = avcn(n[2])
    Avcez = avcn(n[3])

    AvCEx = kron3(Avcez, Avcey, speye(n[1]))
    AvCEy = kron3(Avcez, speye(n[2]), Avcex)
    AvCEz = kron3(speye(n[3]), Avcey, Avcex)

    AveCEMat = [AvCEx; AvCEy; AvCEz]

    return AveCEMat

end

# aveEdge2Cell gets averaging mapping from edge to cell-center
function aveEdge2Cell(n::Vector{Int})

    Avecx = kron3(av(n[3]), av(n[2]), speye(n[1]))
    Avecy = kron3(av(n[3]), speye(n[2]), av(n[1]))
    Avecz = kron3(speye(n[3]), av(n[2]), av(n[1]))

    AveECMat = [Avecx Avecy Avecz]

    return AveECMat

end

# aveNode2Cell gets averaging mapping from node to cell-center
function aveNode2Cell(n::Vector{Int})

    Ancx = av(n[1])
    Ancy = av(n[2])
    Ancz = av(n[3])

    AveNCMat = kron3(Ancz, Ancy, Ancx)

    return AveNCMat

end

# aveCell2Node gets averaging mapping from cell-center to node
function aveCell2Node(n::Vector{Int})

    Avcnx = avcn(n[1])
    Avcny = avcn(n[2])
    Avcnz = avcn(n[3])

    AveCNMat = kron3(Avcnz, Avcny, Avcnx)

    return AveCNMat
end
