#----------------------------------------------------------
#
# Functions for computing mesh geometry
#
#----------------------------------------------------------


# meshGeoFace forms Face area matrix
function meshGeoFace{T<:Real}(xLen::Vector{T},yLen::Vector{T},zLen::Vector{T})

    nx = length(xLen);
    ny = length(yLen);
    nz = length(zLen);

    FaceMat = [diag(kron3(sdiag(zLen),sdiag(yLen),speye(nx+1)));
               diag(kron3(sdiag(zLen),speye(ny+1),sdiag(xLen)));
               diag(kron3(speye(nz+1),sdiag(yLen),sdiag(xLen)))]

    # convert to sparse matrix
    FaceMat = sdiag(FaceMat)

    return FaceMat
end


# meshGeoEdge forms edge length Matrix
function meshGeoEdge{T<:Real}(xLen::Vector{T},yLen::Vector{T},zLen::Vector{T})

    nx = length(xLen);
    ny = length(yLen);
    nz = length(zLen);

    EdgeMat = [diag(kron3(speye(nz+1),speye(ny+1),sdiag(xLen)));
               diag(kron3(speye(nz+1),sdiag(yLen),speye(nx+1)));
               diag(kron3(sdiag(zLen),speye(ny+1),speye(nx+1)))];

    # convert to sparse matrix
    EdgeMat = sdiag(EdgeMat)

    return EdgeMat
end


# meshVolume forms mesh volume matrix
function meshGeoVolume{T<:Real}(xLen::Vector{T},yLen::Vector{T},zLen::Vector{T})

    VolMat = kron3(sdiag(zLen),sdiag(yLen),sdiag(xLen));

    return VolMat

end
