function mxListIntoSparse(MX::Vector{Array{Float64,2}})
    dim = size(MX[1])[1]
    values = [[MXi[i, j] for MXi in MX] for i = 1:dim for j = 1:dim]
    x, y = sparseIndexes(dim, 1)
    return sparse(x, y, values)
end