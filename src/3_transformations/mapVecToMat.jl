function mapVecToMat(V::Array{Float64,1})
    vecLen = length(V)
    n = trunc(Int64, (-1 + sqrt(1 + 8 * vecLen)) / 2)
    Mat = zeros(n, n)
    idx = 0
    for i = 1:n
        Mat[1:n-i+1, i:n] += diagm(V[1+idx:n-i+1+idx])
        idx += n - i + 1
    end
    return LinearAlgebra.Symmetric(Mat)
end