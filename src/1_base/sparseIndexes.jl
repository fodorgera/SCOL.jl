function sparseIndexes(a, n)
    # a: inserted matrix dim
    # n: node num
    xs = zeros((a * n)^2)
    ys = zeros((a * n)^2)

    xrange = vcat([[ai for aii = 1:a] for ai = 1:a]...)
    xrange_row = vcat([xrange for ni = 1:n]...)
    xs = vcat([xrange_row .+ (ni - 1) * a for ni = 1:n]...)
    yrange = vcat([[aii for aii = 1:a] for ai = 1:a]...)
    yrange_row = vcat([yrange .+ (ni - 1) * a for ni = 1:n]...)
    ys = vcat([yrange_row for ni = 1:n]...)
    return xs, ys
end