#----------------------------------------------------------
# difference stencil
#----------------------------------------------------------

# ddx is to form 1D difference sparse matrix from node to center
function ddx(n::Int)
    return spdiagm((-ones(n), ones(n)), (0, 1), n, n+1)
end

# ddxC2N to form 1D difference sparse matrix from grid center to
# node
function ddxC2N(n::Int)

    dC2N = spdiagm((-ones(n), ones(n)), (-1, 0), n+1, n)
    dC2N[end, end] = 1

    return dC2N
end


# av is to form 1D averaging matrix from node to cell-center
function av(n::Int)
    return spdiagm((0.5*ones(n), 0.5*ones(n)), (0, 1), n, n+1)
end

# avcn is to form 1D averaging matrix from cell-center to node
function avcn(n::Int)
    avn = spdiagm((0.5*ones(n), 0.5*ones(n)), (-1, 0), n+1, n)
    avn[1, 1]     = 1
    avn[end, end] = 1

    return avn
end

"""
`avnc(n)` is to form 1D averaging matrix from node to cell-center

"""
function avnc(n::Int)
    avn = spdiagm((0.5*ones(n), 0.5*ones(n)), (0, 1), n, n+1)

    return avn
end

# sdiag is to form sparse diagonal matrix
function sdiag{T<:Number}(v::Vector{T})
    return spdiagm(vec(v), 0, length(v), length(v))
end

# inverse matrix of a sparse diagonal matrix
function sdInv(sdMat::SparseMatrixCSC)

    M = diag(sdMat)
    return sdiag(1./M)

end

# cell-center 1D gradient operator with boundary condition
function ddxCellGradBC(n::Int, BC::String)

    GradBC = spdiagm((-ones(n), ones(n)), (-1,0), n+1, n)
    BC     = lowercase(BC)

    # first side
    if BC == "neumann"
        GradBC[1, 1] = 0
        # second side
        GradBC[end, end]  = 0
    elseif BC == "dirichlet"
        GradBC[1, 1] = -2
        # second side
        GradBC[end, end] = 2
    end

    return GradBC

end

#----------------------------------------------------------
# matrix Utilities
#----------------------------------------------------------

# kronecker tensor products for three vector or matrices
function kron3{T<:SparseMatrixCSC}(A::T, B::T, C::T)

    return kron(A, kron(B, C))
end

#
function ndgrid{T}(xVer::AbstractVector{T}, yVer::AbstractVector{T},
                   zVer::AbstractVector{T})

    nx = length(xVer)
    ny = length(yVer)
    nz = length(zVer)

    xVer = reshape(xVer, nx, 1, 1)
    yVer = reshape(yVer, 1, ny, 1)
    zVer = reshape(zVer, 1, 1, nz)

    idx = ones(Int, nx)
    idy = ones(Int, ny)
    idz = ones(Int, nz)

    return xVer[:, idy, idz], yVer[idx, :, idz], zVer[idx, idy, :]

end
