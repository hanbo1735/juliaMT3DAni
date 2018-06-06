#
# iterSolver solves the linear equations using Krylov methods.
#
function iterSolver(Ke::SparseMatrixCSC, rhs::AbstractArray, PC::Function,
                    linSolver::String, iterLogID::IOStream, tol::Real)


    (nD, ns) = size(rhs)
    xt = zeros(Complex{Float64}, nD, ns)

    # D  = diag(Ke)
    # MM(x) = D.\x
    # MM(x) = tril(Ke)\(D.*(triu(Ke)\x))

    Af(x) = Ke*x

    if uppercase(linSolver) == "BICGSTB"
        #Af(x) = Ke*x
        for j=1:ns
            tic()
            @printf(iterLogID, "%s %3d\n","Pol No.:", j)
            @printf(iterLogID, "%4s\t%7s\n","iter","relres")

            xt[:,j], flag, resLast, iter, resvec = bicgstb(Af, rhs[:,j], tol=tol,
                                                   maxIter=1000, M1=PC, out=2)

            for k=2:length(resvec)
                @printf(iterLogID, "%4d %11.2e\n", k-1, resvec[k])
            end

            elapseTime = toc()

            @printf(iterLogID, "%s %g %s\n","elapsed time:", elapseTime, "seconds.")
        end

    elseif uppercase(linSolver) == "QMR"
        #Af(x) = Ke*x
        for j=1:ns
            tic()
            @printf(iterLogID, "%s %3d\n","Pol No.:", j)
            @printf(iterLogID, "%4s\t%7s\n","iter","relres")

            xt[:,j], flag, resLast, iter, resvec = qmr(Af, rhs[:,j], tol=tol,
                                                   maxIter=1000, M=PC, out=2)

            for k=2:length(resvec)
                @printf(iterLogID, "%4d %11.2e\n", k-1, resvec[k])
            end

            elapseTime = toc()

            @printf(iterLogID, "%s %g %s\n","elapsed time:", elapseTime, "seconds.")
        end

    end # uppercase(linSolver) == ?

    return xt

end

include("qmr.jl")
