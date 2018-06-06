"""
 module `MT3DFileIO` defines composite types and I/O routines for data and model
 parameters for 3D MT problem in general anisotropic media.

 - by Bo Han, Oct 2016

"""
module MT3DFileIO

import YeeGridOperators: CondTensor, EMTensorMesh
import YeeGridOperators: rotateTensor3!, rotateTensor3Inv!
import Utilities.locateNearest

# type
export DataInfo

# functions
export readMT3DData
export writeMT3DData
export readModel3DAni
export writeModel3DAni
export writeMT3DResp
export writeEdgeFields


# This type contains accompanying informations about data, e.g. receiver locations,
# frequencies, data types, but does not contain data (and error) itself.
type DataInfo{T<:Real}

    # receiver locations
    rxLoc::Array{T, 2}

    # frequencies
    freqs::Vector{T}

    # data types (impedance, apparent resistivity, tipper, etc.)
    dataType::AbstractString

    # the specific data components associated with each data type.
    dataComp::Vector{String}

    # receiver index
    rxID::Vector{Int}

    # frequency index
    freqID::Vector{Int}

    # data type index
    dcID::Vector{Int}

end # type


include("readMT3DData.jl")
include("writeMT3DData.jl")
include("readModel3DAni.jl")
include("writeModel3DAni.jl")
include("writeMT3DResp.jl")
include("writeEdgeFields.jl")


end # module
