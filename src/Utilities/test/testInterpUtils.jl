#
# test InterpUtils in module Utilities
# (c) PRH, EOAS, 15 April, 2015
#

using Base.Test
using Utilities

println("=== Testing InterpUtils ===")

xLocI = vec([1 2 3])
yLocI = vec([1 2 3])
zLocI = vec([1 2 3])
xLocF = vec([1.0 2.0 3.0])
yLocF = vec([1.0 2.0 3.0])
zLocF = vec([1.0 2.0 3.0])

p0 = 1.5
p1 = vec([1.5,2.2])
p2 = [1.5 1.5; 2.5 2.5]
p3 = [1.5 1.5 1.5; 2.5 2.5 2.5]

# test function linearInterp
println("Testing function linearInterp ...")
indL, indR, wL, wR = linearInterp(p0, xLocF)
@test (wL, wR) == (0.5, 0.5)

println("Testing function interpMat1D ...")
inds1D, weights1D = interpMat1D(p1, xLocF)

println("Testing function interpMat2D ...")
inds2D, weights2D = interpMat2D(p2, xLocF, yLocF)

println("Testing function InterpMat3D ...")
inds3D,weights3D = interpMat3D(p3, xLocF, yLocF, zLocF)

println("=== all function in InterpUtils passed ===")