#
# test SparseMatUtils in module Utilities
# (c) PRH, EOAS, 15 April, 2015
#

using Base.Test
using Utilities

#
function testSparseUtils()

    println("=== Testing SparseMatUtils ===")

    N = 4
    # test function ddx
    println("Testing function ddx ...")
    DiffN2C = ddx(N);
    @test size(DiffN2C) == (N,N+1)

    # test function ddxC2N
    println("Testing function ddxC2N ...")
    DiffC2N = ddxC2N(N);
    @test size(DiffC2N) == (N+1,N)

    # test function av
    println("Testing function av ...")
    avMat = av(N);

    # test function avcn
    println("Testing function avcn ...")
    avcnMat = avcn(N);

    # test function sdiag
    a = vec([1 2 3 4]);
    b = vec([1.0 2.0 3.0 4.0]);
    println("Testing function sdiag ...")
    aMat = sdiag(a);
    bMat = sdiag(b);

    # test function sdInv
    println("Testing function sdInv ...")
    cMat = sdInv(aMat);
    dMat = sdInv(bMat);

    # test function ddxCellGradBC
    println("Testing function ddxCellGradBC ...")
    neumMat = ddxCellGradBC(N,"neumann")
    diriMat = ddxCellGradBC(N,"dirichlet")

    # test Kron3
    println("Testing function kron3 ...")
    kMat  = kron3(ddx(N),ddx(N),ddx(N))
    k2Mat = kron(ddx(N),kron(ddx(N),ddx(N)))
    @test kMat == k2Mat

    # test ndgrid
    println("Testing function ndgrid ...")


    println("=== All functions in module Utilities passed. ===")

end
