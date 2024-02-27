function buildAddVecFromFunction(VEC::Function, j, τ, Coll)
    nodes = Coll.nodes
    div = Coll.div
    order = size(VEC(0.0))[1]
    xs = Coll.xs
    ls = Coll.ls
    mult = Coll.mult
    subdiv = innerSubdiv
    values_k = [ItoSVec(0, SVector{subdiv + 1}([0.0 for i = 1:subdiv+1])) for j = 1:(nodes*order)^2]
    values = [deepcopy(values_k) for k = 1:div]
    for k = 1:div
        for i = 1:nodes
            if i > 1
                # dη=abs(xs[i-1]-xs[i])/subdiv;
                # xvals=xs[i-1]:dη:xs[i];
                xvals = range(xs[i-1], xs[i], length = subdiv + 1)
                CoeffMX = mxListIntoSparse(diagm.(VEC.(tjk(j, k, τ, Coll).(xvals))))
                for o = 1:order
                    for p = 1:order
                        values[k][(i-1)*nodes*order^2+(o-1)*order*nodes+(i-1)*order+p].v = SVector{subdiv + 1}(sqrt(τ * mult / div * abs(xs[i] - xs[i-1])) .* CoeffMX[o, p])
                        values[k][(i-1)*nodes*order^2+(o-1)*order*nodes+(i-1)*order+p].idx = parse(Int64, "$(j)$(k)$(i)")
                    end
                end
            end
        end
    end

    values_full = vcat(values...)
    return values_full
end