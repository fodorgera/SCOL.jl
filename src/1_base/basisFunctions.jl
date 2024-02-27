function legendreFunctions(_nodes::Int64)
    nodes = _nodes - 1
    over(n::Int64, k::Int64) = factorial(n) / (factorial(k) * factorial(n - k))
    xs = -1.0:2/(nodes):1.0
    Le = [x -> sum([over(n, k) * over(n + k, k) * ((x - 1) / 2)^k for k = 0:n]) for n = 0:nodes]
    return xs, Le
end

function lagrangeFunctions(nodes::Int64)
    xs = 0.0:1/(nodes-1):1.0
    ls = [x -> prod([(x - xs[i]) / (xs[j] - xs[i]) for i in 1:length(xs) if i != j]) for j = 1:nodes]
    return xs, ls
end

function chebyshevFunctions(_nodes::Int64)
    nodes = _nodes - 1
    Ch = [x -> cos(i * acos(x)) for i = 0:nodes]
    cosxs = [cos((nodes - j) * π / nodes) for j = 0:nodes]
    return cosxs, Ch
end

function lagrangeChebyshevFunctions(_nodes::Int64)
    nodes = _nodes - 1
    xs = [cos((nodes - j) * π / nodes) for j = 0:nodes]
    ls = [x -> prod([(x - xs[i]) / (xs[j] - xs[i]) for i in 1:length(xs) if i != j]) for j = 1:_nodes]
    return xs, ls
end

function selectMethod(method::String, nodes::Int64)
    if method == "Lagrange"
        return 1.0, lagrangeFunctions(nodes)
    elseif method == "Legendre"
        return 1 / 2, legendreFunctions(nodes)
    elseif method == "Chebyshev"
        return 1 / 2, chebyshevFunctions(nodes)
    elseif method == "LagrangeChebyshev" || method == "ChebyshevLagrange"
        return 1 / 2, lagrangeChebyshevFunctions(nodes)
    else
        println("unknown method")
    end
end