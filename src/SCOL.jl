module SCOL

using QuadGK, LinearAlgebra, SemiDiscretizationMethod, MDBM, StaticArrays, InteractiveUtils, Arpack

export getMappingOfSystem

include("./dependencies.jl")
include("./1_base/index.jl")
include("./2_assembly/index.jl")
include("./3_transformations/index.jl")
include("./4_build/index.jl")
include("./5_periodic/index.jl")
include("./6_moments/index.jl")
include("./7_evaluation/index.jl")

end
