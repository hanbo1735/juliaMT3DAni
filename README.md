# juliaMT3DAni
`juliaMT3DAni` is a 3D forward modeling program for MT surveys in general anisotropic media written in [the Julia language](https://julialang.org/).

For the details regarding the algorithm and implementation, please refer to: 
> Han, B., Y. Li, and G. Li, 2018, 3D forward modeling of magnetotelluric fields in
> general anisotropic media and its numerical implementation in Julia:
> Geophysics, 83(4), F29-F40; DOI:
> [10.1190/geo2017-0515.1](http://doi.org/10.1190/geo2017-0515.1).


## Prerequisite
`juliaMT3DAni` utilizes three third-party Julia packages:
[KrylovMethods.jl](https://github.com/lruthotto/KrylovMethods.jl), [MUMPS.jl](https://github.com/JuliaSparse/MUMPS.jl) and [Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl).
