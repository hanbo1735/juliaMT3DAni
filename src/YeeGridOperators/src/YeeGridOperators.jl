"""
 module `YeeGridOpterators` defines operators for Yee's tensor product grid.

 - Originally written by Ronghua Peng, University of British Columbia,
   22 April, 2015; modified by Bo Han, Ocean University of China, Oct 2016.

"""

module YeeGridOperators

using Utilities

# you must import them explicitly if you want to extend them
import Utilities: rotateTensor3!, rotateTensor3Inv!

export setupOperators!
export CondTensor, EMTensorMesh
export rotateTensor3!, rotateTensor3Inv!

type CondTensor{T<:Real}
  sig::Array{T,2}
  rotated::Bool
end # TensorSigma

type EMTensorMesh{T<:Real}

    # width of cells in x direction
    xLen::Vector{T}
    # width of cells in y direction
    yLen::Vector{T}
    # width of cells in z direction excluding airlayer
    zLen::Vector{T}
    # airlayer
    airLayer::Vector{T}
    # cell numbers in x,y,z direction
    gridSize::Vector{Int}
    # origin of mesh
    origin::Vector{T}
    # total number of cells
    nGrid::Int
    # conductivity of mesh
    sigma::CondTensor
    # Volume of cells
    Vol::SparseMatrixCSC
    # edge curl
    Curl::SparseMatrixCSC
    # nodal gradient
    Grad::SparseMatrixCSC
    # face divergence
    Div::SparseMatrixCSC
    # averaging mapping from edge to cell-center
    AveEC::SparseMatrixCSC
    # averaging mapping from face to cell-center
    AveFC::SparseMatrixCSC
    # averaging mapping from node to cell-center
    AveNC::SparseMatrixCSC

end # EMTensorMesh

#
# rotateTensor3!
#
function rotateTensor3!(sigma::CondTensor)
    rotateTensor3!(sigma.sig)
    sigma.rotated = true
    return sigma
end

#
# rotateTensor3Inv!
#
function rotateTensor3Inv!(sigma::CondTensor)
    rotateTensor3Inv!(sigma.sig)
    sigma.rotated = false
    return sigma
end

#
function setupOperators!(emMesh::EMTensorMesh)
    xLen = emMesh.xLen
    yLen = emMesh.yLen
    zLen = emMesh.zLen
    gridSize = emMesh.gridSize

    emMesh.Curl  = getEdgeCurl(xLen, yLen, zLen)
    emMesh.Div   = getFaceDivergence(xLen, yLen, zLen)
    emMesh.Grad  = getNodalGradient(xLen, yLen, zLen)
    emMesh.AveEC = aveEdge2Cell(gridSize)
    emMesh.AveFC = aveFace2Cell(gridSize)
    emMesh.AveNC = aveNode2Cell(gridSize)
    emMesh.Vol   = meshGeoVolume(xLen, yLen, zLen)
end


include("linearOperators.jl")
include("meshIO.jl")

end # module
