"""
`getAnisoOperator` gets the coefficients arising from off-diagonal anisotropy.

Input:
    gridSize  :: Vector            - cell numbers in three directions.
    Vol       :: SparseMatrixCSC   - cell volume matrix.
    sigxy, sigxz, sigyz :: Vector  - off-diagonal entries of the conductivity tensor.

Output:
    A :: SparseMatrixCSC   - the coefficients arising from off-diagonal anisotropy.

"""

function getAnisoOperator{T<:Real}(gridSize::Vector{Int}, Vol::SparseMatrixCSC,
                                   sigxy::Vector{T}, sigxz::Vector{T}, sigyz::Vector{T})


    nx = gridSize[1]; ny = gridSize[2]; nz = gridSize[3]
    nEx = nx*(ny+1)*(nz+1)
    nEy = (nx+1)*ny*(nz+1)
    nEz = (nx+1)*(ny+1)*nz

    # define some utility anonymous functions
    av   = (n::Int) -> spdiagm((0.5*ones(n), 0.5*ones(n)), (0, 1), n, n+1)
    oeye = (n::Int) -> vcat(spzeros(1,n), speye(n))     # size: (n+1,n)
    eyeo = (n::Int) -> vcat(speye(n), spzeros(1,n))

    # Ex -> Ey
    a_xy1 = kron(speye(nz+1), kron(oeye(ny), av(nx)))  # size: (nEx,nEy)
    a_xy2 = kron(speye(nz+1), kron(eyeo(ny), av(nx)))
    A_xy1 = [spzeros(nEx,nEx)    a_xy1    spzeros(nEx,nEz)]
    A_xy2 = [spzeros(nEx,nEx)    a_xy2    spzeros(nEx,nEz)]

    c2E_xy1 = kron((2*av(nz))', kron(oeye(ny), speye(nx)))  # size: (nEx,nCell)
    c2E_xy2 = kron((2*av(nz))', kron(eyeo(ny), speye(nx)))
    Msigxy1 = spdiagm(c2E_xy1 * 0.25 * Vol * sigxy)
    Msigxy2 = spdiagm(c2E_xy2 * 0.25 * Vol * sigxy)

    Axy = Msigxy1 * A_xy1 + Msigxy2 * A_xy2  # size: (nEx,nE)


    # Ex -> Ez
    a_xz1 = kron(oeye(nz), kron(speye(ny+1), av(nx)))  # size: (nEx,nEz)
    a_xz2 = kron(eyeo(nz), kron(speye(ny+1), av(nx)))
    A_xz1 = [spzeros(nEx,nEx)    spzeros(nEx,nEy)    a_xz1]
    A_xz2 = [spzeros(nEx,nEx)    spzeros(nEx,nEy)    a_xz2]

    c2E_xz1 = kron(oeye(nz), kron((2*av(ny))', speye(nx))) # size: (nEx,nCel)
    c2E_xz2 = kron(eyeo(nz), kron((2*av(ny))', speye(nx)))
    Msigxz1 = spdiagm(c2E_xz1 * 0.25 * Vol * sigxz)
    Msigxz2 = spdiagm(c2E_xz2 * 0.25 * Vol * sigxz)

    Axz = Msigxz1 * A_xz1 + Msigxz2 * A_xz2  # size: (nEx,nE)

    AX = Axy + Axz


    # Ey -> Ex
    a_yx1 = kron(speye(nz+1), kron(av(ny), oeye(nx)))  # size: (nEy,nEx)
    a_yx2 = kron(speye(nz+1), kron(av(ny), eyeo(nx)))
    A_yx1 = [a_yx1    spzeros(nEy,nEy)    spzeros(nEy,nEz)]
    A_yx2 = [a_yx2    spzeros(nEy,nEy)    spzeros(nEy,nEz)]

    c2E_yx1 = kron((2*av(nz))', kron(speye(ny), oeye(nx)))  # size: (nEy,nCell)
    c2E_yx2 = kron((2*av(nz))', kron(speye(ny), eyeo(nx)))
    Msigyx1 = spdiagm(c2E_yx1 * 0.25 * Vol * sigxy)
    Msigyx2 = spdiagm(c2E_yx2 * 0.25 * Vol * sigxy)

    Ayx = Msigyx1 * A_yx1 + Msigyx2 * A_yx2  # size: (nEy,nE)


    # Ey -> Ez
    a_yz1 = kron(oeye(nz), kron(av(ny), speye(nx+1)))  # size: (nEy,nEz)
    a_yz2 = kron(eyeo(nz), kron(av(ny), speye(nx+1)))
    A_yz1 = [spzeros(nEy,nEx)    spzeros(nEy,nEy)    a_yz1]
    A_yz2 = [spzeros(nEy,nEx)    spzeros(nEy,nEy)    a_yz2]

    c2E_yz1 = kron(oeye(nz), kron(speye(ny), (2*av(nx))')) # size: (nEy,nCell)
    c2E_yz2 = kron(eyeo(nz), kron(speye(ny), (2*av(nx))'))
    Msigyz1 = spdiagm(c2E_yz1 * 0.25 * Vol * sigyz)
    Msigyz2 = spdiagm(c2E_yz2 * 0.25 * Vol * sigyz)

    Ayz = Msigyz1 * A_yz1 + Msigyz2 * A_yz2 # size: (nEx,nE)

    AY = Ayx + Ayz


    # Ez -> Ex
    a_zx1 = kron(av(nz), kron(speye(ny+1), oeye(nx)))  # size: (nEz,nEx)
    a_zx2 = kron(av(nz), kron(speye(ny+1), eyeo(nx)))
    A_zx1 = [a_zx1    spzeros(nEz,nEy)    spzeros(nEz,nEz)]
    A_zx2 = [a_zx2    spzeros(nEz,nEy)    spzeros(nEz,nEz)]

    c2E_zx1 = kron(speye(nz), kron((2*av(ny))', oeye(nx)))  # size: (nEz,nCell)
    c2E_zx2 = kron(speye(nz), kron((2*av(ny))', eyeo(nx)))
    Msigzx1 = spdiagm(c2E_zx1 * 0.25 * Vol * sigxz)
    Msigzx2 = spdiagm(c2E_zx2 * 0.25 * Vol * sigxz)

    Azx = Msigzx1 * A_zx1 + Msigzx2 * A_zx2  # size: (nEz,nE)


    # Ez -> Ey
    a_zy1 = kron(av(nz), kron(oeye(ny), speye(nx+1)))  # size: (nEz,nEy)
    a_zy2 = kron(av(nz), kron(eyeo(ny), speye(nx+1)))
    A_zy1 = [spzeros(nEz,nEx)    a_zy1    spzeros(nEz,nEz)]
    A_zy2 = [spzeros(nEz,nEx)    a_zy2    spzeros(nEz,nEz)]

    c2E_zy1 = kron(speye(nz), kron(oeye(ny), (2*av(nx))'))  # size: (nEz,nCell)
    c2E_zy2 = kron(speye(nz), kron(eyeo(ny), (2*av(nx))'))
    Msigzy1 = spdiagm(c2E_zy1 * 0.25 * Vol * sigyz)
    Msigzy2 = spdiagm(c2E_zy2 * 0.25 * Vol * sigyz)

    Azy = Msigzy1 * A_zy1 + Msigzy2 * A_zy2  # size: (nEz,nE)

    AZ = Azx + Azy

    return vcat(AX, AY, AZ)

end
