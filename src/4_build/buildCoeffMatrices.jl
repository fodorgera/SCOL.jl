function buildCoeffMatrices(A::Array{Float64,2}, B::Array{Float64,2}, c::Array{Float64,1}, α::Array{Float64,2}, β::Array{Float64,2}, σ::Array{Float64,1}, τ::Float64, collMethod::CollocationMethod{xsT,lsT}) where {xsT,lsT}
    F0, F1, cMat = assembleDeterministic(A, B, c, τ, collMethod)
    G0, G1, σMat = assembleStochastic(α, β, σ, τ, collMethod)
    return CoeffMatrices(F0, F1, cMat, G0, G1, σMat)
end