function mxmult(a::Array{Float64,2}, b::Array{Float64,2})
    sa = size(a)[1]
    sb = size(b)[1]
    mx = zeros(sa * sb, sa * sb)
    for i in 1:sb
        for j in 1:sb
            mx[(i-1)*sa+1:i*sa, (j-1)*sa+1:j*sa] = a .* b[i, j]
        end
    end
    return mx
end

function mxmult(a::Array{Float64,2}, b::SparseMatrixCSC{ItoSVec{SArray{Tuple{innerSubdiv + 1},Float64,1,innerSubdiv + 1},Float64},Int64})
    sa = size(a)[1]
    sb = size(b)[1]
    mx = spzeros(ItoSVec{SArray{Tuple{innerSubdiv + 1},Float64,1,innerSubdiv + 1},Float64}, sa * sb, sa * sb)
    for i in 1:sb
        for j in 1:sb
            for k in 1:sa
                for l in 1:sa
                    try
                        mx[(i-1)*sa+k, (j-1)*sa+l] = ItoSVec(b[i, j].idx, a[k, l] * b[i, j].v)
                    catch e
                    end
                end
            end
        end
    end
    return mx
end

function mxmult(A::Array{Float64,2}, values::Array{ItoSVec{SArray{Tuple{innerSubdiv + 1},Float64,1,innerSubdiv + 1},Float64},1})
    # B is a sparse matrix
    nA = length(A[1, :])
    nval = length(values)
    subdiv = innerSubdiv
    n = nA^2 * nval
    new_values = [ItoSVec(0, SVector{subdiv + 1}([0.0 for i = 1:subdiv+1])) for z = 1:n]
    for i = 1:nval
        num = 0
        for row = 1:nA
            for col = 1:nA
                num += 1
                new_values[nA^2*(i-1)+num].v = A[row, col] .* values[i].v
                new_values[nA^2*(i-1)+num].idx = values[i].idx
            end
        end
    end

    return new_values
end