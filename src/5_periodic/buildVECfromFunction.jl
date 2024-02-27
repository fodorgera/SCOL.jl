function buildVECfromFunction(VEC::Function, j, τ, Coll)
    order = size(VEC(0.0))[1]
    VEC_k = zeros(Coll.nodes * order)
    VEC_j = [deepcopy(VEC_k) for k = 1:Coll.div]
    dim = Coll.nodes * Coll.div * order
    for k = 1:Coll.div
        VEC_j[k][order+1:end] = hcat([quadgk(x -> τ * Coll.mult / Coll.div * VEC(tjk(j, k, τ, Coll)(x)), Coll.xs[i], Coll.xs[i+1], atol = 0.0001)[1] for i = 1:Coll.nodes-1]...)
    end
    VEC = zeros(dim)
    for k = 1:Coll.div
        VEC[(k-1)*Coll.nodes*order+1:k*Coll.nodes*order] = VEC_j[k]
    end
    return VEC
end