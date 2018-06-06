"""
 module `LinearSolver` defines routines to solve large sparse linear systems.

 - by Bo Han, Oct 2016

"""
module LinearSolver

using MUMPS
using Pardiso
using KrylovMethods

export mumpsSolver, pardiSolver, iterSolver
export qmr

include("mumpsSolver.jl")
include("iterSolver.jl")
include("pardiSolver.jl")

end # LinearSolver
