function buildSparseFromFunction2(MX::Function, j, τ, Coll)
    nodes = Coll.nodes
    div = Coll.div
    order = size(MX(0.0))[1]
    xs = Coll.xs
    ls = Coll.ls
    mult = Coll.mult

    subdiv = innerSubdiv
    values_k = [ItoSVec(0, SVector{subdiv + 1}([0.0 for i = 1:subdiv+1])) for j = 1:(nodes)^2]
    mult_values_k = [ItoSVec(0, SVector{subdiv + 1}([0.0 for i = 1:subdiv+1])) for z = 1:(nodes*order)^2]
    values = [deepcopy(values_k) for k = 1:div]
    mult_values = [deepcopy(mult_values_k) for k = 1:div]

    # dη=abs(xs[2]-xs[1])/subdiv;
    # xvals=xs[1]:dη:xs[2];
    xvals = range(xs[1], xs[2], length = subdiv + 1)
    TemplateCoeffMX = mxListIntoSparse(MX.(tjk(j, 1, τ, Coll).(xvals)) * 0.0)
    coeffMXs_k = [deepcopy(TemplateCoeffMX) for i = 1:nodes^2]
    coeffMXs = [deepcopy(coeffMXs_k) for k = 1:div]
    for k = 1:div
        for i = 1:nodes
            if i > 1
                # dη=abs(xs[i-1]-xs[i])/subdiv;
                # xvals=xs[i-1]:dη:xs[i];
                xvals = range(xs[i-1], xs[i], length = subdiv + 1)
                CoeffMX = mxListIntoSparse(MX.(tjk(j, k, τ, Coll).(xvals)))
                for l = 1:nodes
                    CoeffMXVal = [CoeffMX[s, z] for s = 1:order for z = 1:order]
                    coeffMXs[k][nodes*(i-1)+l] = CoeffMX
                    values[k][nodes*(i-1)+l].v = SVector{subdiv + 1}(sqrt(τ * mult / div * abs(xs[i] - xs[i-1])) .* ls[l].(xvals))
                    values[k][nodes*(i-1)+l].idx = parse(Int64, "$(j)$(k)$(i)")
                end
            end
        end


        for z = 1:nodes^2
            mx = coeffMXs[k][z]
            # println(mx)
            num = 0
            for row = 1:order
                for col = 1:order
                    num += 1
                    mult_values[k][order^2*(z-1)+num].v = mx[row, col] .* values[k][z].v
                    mult_values[k][order^2*(z-1)+num].idx = values[k][z].idx
                end
            end
        end
    end

    values_full = vcat(mult_values...)
    return values_full
end