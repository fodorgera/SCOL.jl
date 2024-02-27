struct DetMXs
    F0::Array{Float64,2}
    F1::Array{Float64,2}
    ct::Array{Float64,2}
end

struct StochMXs
    αt::SparseMatrixCSC{ItoSVec{StaticArrays.SArray{Tuple{innerSubdiv + 1},Float64,1,innerSubdiv + 1},Float64},Int64}
    βt::SparseMatrixCSC{ItoSVec{StaticArrays.SArray{Tuple{innerSubdiv + 1},Float64,1,innerSubdiv + 1},Float64},Int64}
    σt::SparseMatrixCSC{ItoSVec{StaticArrays.SArray{Tuple{innerSubdiv + 1},Float64,1,innerSubdiv + 1},Float64},Int64}
end