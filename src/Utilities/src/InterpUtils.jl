#
# InterpUtils defines some routines for linear interpolation
# (c) PRH, EOAS, 15 April, 2015
#
#  Modified by Bo Han, July 2016
#  rename some functions, and add some functions.

#
# LinearInterp gets 1D linear interpolation points and weights for a single point
# 
function linearInterp{T<:Real}(point::T, x::Vector{T})

    val, ind = findmin(abs.(point - x))
    # location of point
    if point - x[ind] > 0
    # point on the right
        indL = ind
        indR = ind + 1
    else
    # point on the left
        indL = ind - 1
        indR = ind
    end

    # ensure interpolation points within the bound
    n = length(x)
    indL = maximum([minimum([indL, n]), 1])
    indR = maximum([minimum([indR, n]), 1])

    if indL == indR
        return indL, indR, 0.5, 0.5
    end

    # interpolation weights
    xLen = x[indR] - x[indL]
    wL   = 1 - (point - x[indL]) / xLen
    wR   = 1 - (x[indR] - point) / xLen

    return indL, indR, wL, wR

end


#
# linearInterpMat gets 1D linear interpolation points and weights
# and use them to construct the interpolation matrix.
#
function linearInterpMat{T<:Real}(point::Vector{T}, x::Vector{T}, getMat=true)

    npts  = length(point)
    nNode = length(x)

    inds    = cell(npts)
    weights = zeros(Float64, npts, 2)
    interpMat = spzeros(Float64, nNode, npts)

    for i = 1:npts
        indL, indR, wL, wR = linearInterp(point[i], x)
        inds[i] = (indL,indR)
        weights[i,:] = [wL wR]
        interpMat[:,i] = sparsevec([indL;indR], [wL;wR], nNode)
    end

    if !getMat
        return inds, weights
    end

    return interpMat

end


#
# bilinearInterpMat gets bilinear interpolation points and weights
# and use them to construct the interpolation matrix.
#
function bilinearInterpMat{T<:Real}(point::Array{T,2}, x::Vector{T},
                                    y::Vector{T}, getMat=true)

    npts = size(point, 1)
    nx = length(x);  ny = length(y)
    nNode = nx*ny

    inds    = Array{Any}(npts)
    weights = zeros(Float64, npts, 4)
    interpMat = spzeros(Float64, nNode, npts)

    for i = 1:npts
        indxL, indxR, wxL, wxR = linearInterp(point[i, 1], x)
        indyL, indyR, wyL, wyR = linearInterp(point[i, 2], y)
        inds[i] = [(indxL, indyL), (indxR, indyL),
                   (indxL, indyR), (indxR, indyR)]

        w1 = wxL * wyL
        w2 = wxR * wyL
        w3 = wxL * wyR
        w4 = wxR * wyR

        w = [w1 w2 w3 w4]
        weights[i,:] = w

        idx = zeros(Int, 4)
        for j = 1:4
            idx[j] = nx * (inds[i][j][2]-1) + inds[i][j][1]
        end
        interpMat[:,i] = sparsevec(idx, vec(w), nNode)

    end

    if !getMat
        return inds, weights
    end

    return interpMat

end


#
# trilinearInterpMat gets trilinear interpolation points and weights
# and use them to construct the interpolation matrix.
#
function trilinearInterpMat{T<:Real}(point::Array{T, 2}, x::Vector{T},
                                     y::Vector{T}, z::Vector{T}, getMat=true)

    npts = size(point,1)
    nx = length(x);  ny = length(y);  nz = length(z)
    nNode = nx*ny*nz

    inds    = cell(npts)
    weights = zeros(Float64, npts, 8)
    interpMat = spzeros(Float64, nNode, npts)

    for i = 1:npts
        indxL, indxR, wxL, wxR = linearInterp(point[i, 1], x)
        indyL, indyR, wyL, wyR = linearInterp(point[i, 2], y)
        indzL, indzR, wzL, wzR = linearInterp(point[i, 3], z)
        inds[i] = [(indxL, indyL, indzL),
                   (indxR, indyL, indzL),
                   (indxL, indyR, indzL),
                   (indxR, indyR, indzL),
                   (indxL, indyL, indzR),
                   (indxR, indyL, indzR),
                   (indxL, indyR, indzR),
                   (indxR, indyR, indzR)]

        w1 = wxL * wyL * wzL
        w2 = wxR * wyL * wzL
        w3 = wxL * wyR * wzL
        w4 = wxR * wyR * wzL
        w5 = wxL * wyL * wzR
        w6 = wxR * wyL * wzR
        w7 = wxL * wyR * wzR
        w8 = wxR * wyR * wzR

        w = [w1 w2 w3 w4 w5 w6 w7 w8]
        weights[i,:] = w

        idx = zeros(Int, 8)
        for j = 1:8
            idx[j] = nx * ny * (inds[i][j][3]-1) + nx * (inds[i][j][2]-1) + inds[i][j][1]
        end
        interpMat[:,i] = sparsevec(idx, vec(w), nNode)

    end

    if !getMat
        return inds, weights
    end

    return interpMat

end


#
# locateNearest finds the location of points to the nearest coordinate location
# (node,cell-center, etc.)
#
function locateNearest{T<:Real}(point::Array{T,2}, xLoc::Vector{T},
                                yLoc::Vector{T}, zLoc::Vector{T})

    location = ones(Int,3)

    # x location
    val,ind = findmin(abs(point[1] - xLoc))
    location[1] = ind

    # y location
    val,ind = findmin(abs(point[2] - yLoc))
    location[2] = ind

    # z location
    val,ind = findmin(abs(point[3] - zLoc))
    location[3] = ind

    return location

end


"""
`locateSegment1D` determines the segment a point reside in.
Note that if the point is at the boundary, assume it's in the left segment.
"""
function locateSegment1D{T<:Real}(xp::T, x::Vector{T})

    n = length(x)
    val, ind = findmin(abs(x - xp))
    if xp <= x[ind]; ind -= 1; end
    ind = maximum([minimum([ind, n-1]), 1])

end
