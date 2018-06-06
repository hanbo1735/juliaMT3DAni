#
# mumpsSolver is the interface to MUMPS.jl.
#
function mumpsSolver(Ke::SparseMatrixCSC, rhs::AbstractArray; ooc::Int=0, saveFac::Bool=false)

    # the coefficient matrix is simply complex symmetric (not hermitan), so an LDLᴴ
    # decomposition can be used. MUMPS package is used for this purpose.

    sym  = 1
    @time Ainv = factorMUMPS(Ke, sym, ooc)

    nD = size(rhs, 1)
    ns = size(rhs, 2)
    xt  = zeros(Complex{Float64}, nD, ns)

    if maximum(abs.(rhs)) == 0.0
        println("All source terms are zeros.")
    else
        xt = applyMUMPS(Ainv, rhs)
    end

    saveFac && ( return xt, Ainv )

    destroyMUMPS(Ainv)
    return xt

end
