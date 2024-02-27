function mapMatToVec(M::Array{Float64,2})
    sizeM = size(M)[1]
    diags = vcat([diag(M[1:end-i+1, i:end]) for i = 1:sizeM]...)
    return diags
end