using KrylovMethods
using Base.Test
using LinearSolver


println("Testing QMR for complex symmetric matrix")
A1 = sprandn(100,100,.1)
A1 = A1 * A1'
A  = A1 + 10*speye(100) + im*(10*speye(100) )
D  = diag(A)
Af = x -> A*x
PC = x -> D.\x
rhs = randn(100) + 1im * randn(100)

x2 = qmr(Af,rhs,tol=1e-6,out=2)
x3 = qmr(Af,rhs,tol=1e-6,maxIter=200,x0=randn(size(rhs))+im*randn(size(rhs)),out=2)
x4 = qmr(Af,rhs,tol=1e-6,maxIter=200,M=PC,out=2)

@test norm(A*x2[1]-rhs)/norm(rhs) < 1e-6
@test norm(A*x3[1]-rhs)/norm(rhs) < 1e-6
@test norm(A*x4[1]-rhs)/norm(rhs) < 1e-6

println("=== QMR: All tests passed ===")
