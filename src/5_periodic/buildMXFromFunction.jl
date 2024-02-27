function buildMXfromFunction(MX::Function, j, τ, Coll)
    order = size(MX(0.0))[1]
    MX_k = zeros(Coll.nodes * order, Coll.nodes * order)
    MX_j = [deepcopy(MX_k) for k = 1:Coll.div]
    dim = Coll.nodes * Coll.div * order
    for k = 1:Coll.div
        MX_j[k][order+1:end, :] = vcat([hcat([quadgk(x -> τ * Coll.mult / Coll.div * MX(tjk(j, k, τ, Coll)(x)) * Coll.ls[l](x), Coll.xs[i], Coll.xs[i+1], atol = 0.0001)[1] for l = 1:Coll.nodes]...) for i = 1:Coll.nodes-1]...)
    end
    MX = zeros(dim, dim)
    for k = 1:Coll.div
        MX[(k-1)*Coll.nodes*order+1:k*Coll.nodes*order, (k-1)*Coll.nodes*order+1:k*Coll.nodes*order] = MX_j[k]
    end
    return MX
end