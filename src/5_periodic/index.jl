include("structs.jl")
include("tjk.jl")

# Deterministic
include("buildMXFromFunction.jl")
include("buildVECfromFunction.jl")
include("buildDetMXsWholeState.jl")

# Stochastic
include("mxListIntoSparse.jl")
include("buildAddVecFromFunction.jl")
include("buildSparseFromFunction2.jl")
include("buildStochMXsWholeState.jl")
