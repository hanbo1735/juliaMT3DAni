#
# pardiSolver is the interface to Pardiso.jl.
#
function pardiSolver(Ke::SparseMatrixCSC, rhs::AbstractArray)

    ps = MKLPardisoSolver()

    set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
    set_nprocs!(ps, 1)

    # you can solve it in a single call
    #@time xt = solve(ps, Ke, rhs, :N)

    # it would be better to do it step by step
    set_matrixtype!(ps, Pardiso.COMPLEX_SYM)
    pardisoinit(ps)

    # out-of-core
    # set_iparm!(ps, 3, 1)
    # set_iparm!(ps, 60, 2)

    A_pardiso = get_matrix(ps, Ke, :N)

    # Analyze the matrix and compute a symbolic factorization &
    # compute the numeric factorization.
    set_phase!(ps, Pardiso.ANALYSIS_NUM_FACT)
    @time pardiso(ps, A_pardiso, rhs)

    # Compute the solutions X using the symbolic factorization.
    set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)
    xt = similar(rhs)
    pardiso(ps, xt, A_pardiso, rhs)

    # Free the PARDISO data structures.
    set_phase!(ps, Pardiso.RELEASE_ALL)
    pardiso(ps)

    return xt

end
