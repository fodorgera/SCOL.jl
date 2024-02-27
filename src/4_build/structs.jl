struct CoeffMatrices
    F0::Array{Float64,2}
    F1::Array{Float64,2}
    cMat::Array{Float64,2}
    G0::SparseMatrixCSC{ItoSVec{SArray{Tuple{innerSubdiv + 1},Float64,1,innerSubdiv + 1},Float64},Int64}
    G1::SparseMatrixCSC{ItoSVec{SArray{Tuple{innerSubdiv + 1},Float64,1,innerSubdiv + 1},Float64},Int64}
    σMat::SparseMatrixCSC{ItoSVec{SArray{Tuple{innerSubdiv + 1},Float64,1,innerSubdiv + 1},Float64},Int64}
end
