"""
 module `Utilities` defines some utility routines for sparse matrix
 manipulation, linear interpolation, etc.

 - by Ronghua Peng, University of British Columbia, 13 April, 2015;
   modified by Bo Han, Ocean University of China, Oct 2016.

"""
module Utilities

# SparseMatUtils
export ddx,ddxC2N
export av, avcn, avnc
export sdiag, sdInv
export ddxCellGradBC
export kron3, ndgrid

# InterpUtils
export linearInterp
export linearInterpMat
export bilinearInterpMat
export trilinearInterpMat
export locateSegment1D
export locateNearest

# TensorUtils
export rotateTensor3!, rotateTensor3Inv!

include("SparseMatUtils.jl")
include("InterpUtils.jl")
include("TensorUtils.jl")

end # Utilities
