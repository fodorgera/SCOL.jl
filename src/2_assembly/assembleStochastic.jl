function assembleStochastic(α::Array{Float64,2}, β::Array{Float64,2}, σ::Array{Float64,1}, τ::Float64, method::CollocationMethod{xsT,lsT}) where {xsT,lsT}
    xs = method.xs
    ls = method.ls
    mult = method.mult
    div = method.div
    nodes = method.nodes

    order = length(α[1, :])
    subdiv = innerSubdiv
    values = [ItoSVec(0, SVector{subdiv + 1}([0.0 for i = 1:subdiv+1])) for j = 1:nodes^2]
    valuesAdd = [ItoSVec(0, SVector{subdiv + 1}([0.0 for i = 1:subdiv+1])) for j = 1:nodes^2]
    for i = 1:nodes
        if i != 1
            # dη=abs(xs[i-1]-xs[i])/subdiv;
            # xvals=xs[i-1]:dη:xs[i];
            xvals = range(xs[i-1], xs[i], length=subdiv + 1)
            valuesAdd[nodes*(i-1)+i].v = SVector{subdiv + 1}([sqrt(τ * mult / div * abs(xs[i] - xs[i-1])) for xi in xvals])
            valuesAdd[nodes*(i-1)+i].idx = i
            for j = 1:nodes
                # dη=abs(xs[i-1]-xs[i])/subdiv;
                # xvals=xs[i-1]:dη:xs[i];
                xvals = range(xs[i-1], xs[i], length=subdiv + 1)
                values[nodes*(i-1)+j].v = SVector{subdiv + 1}(sqrt(τ * mult / div * abs(xs[i] - xs[i-1])) .* ls[j].(xvals))
                values[nodes*(i-1)+j].idx = i
            end
        end
    end

    # the matrix containing the stoch integral vectors
    # getting the indexes
    x_ind, y_ind = sparseIndexes(order, nodes) # order: the dim of the α matrix, nodes: number of trial functions used
    G0_values = mxmult(α, values)
    G1_values = mxmult(β, values)
    σmat = diagm(σ)
    GAdd_values = mxmult(σmat, valuesAdd)
    G0 = sparse(x_ind, y_ind, G0_values)
    G1 = sparse(x_ind, y_ind, G1_values)
    GAdd = sparse(x_ind, y_ind, GAdd_values)
    # calculating according to the divs

    x_ind_A = deepcopy(x_ind)
    x_ind_B = deepcopy(x_ind) .+ nodes * order * (div - 1)
    y_ind_A = deepcopy(y_ind)
    y_ind_B = deepcopy(y_ind) .+ nodes * order * (div - 1)

    G0_nulls = [ItoSVec(0, SVector{subdiv + 1}([0.0 for i = 1:subdiv+1])) for j = 1:(nodes*order)^2]
    # for i in 1:nodes*order
    #     for j in 1:nodes*order
    #         G0_values[nodes*order*(i-1)+j].v = G0[i,j].v;
    #         G0_values[nodes*order*(i-1)+j].idx = G0[i,j].idx;
    #     end
    # end
    #  G0_values = vcat([[ItoSVec(G0[i,j].idx,G0[i,j].v) for j in 1:nodes*order] for i in 1:nodes*order]...);
    #  G1_values = vcat([[ItoSVec(G1[i,j].idx,G1[i,j].v) for j in 1:nodes*order] for i in 1:nodes*order]...);
    #  GAdd_values = vcat([[ItoSVec(GAdd[i,j].idx,GAdd[i,j].v) for j in 1:nodes*order] for i in 1:nodes*order]...);

    function G_gen(div)
        if div > 1
            G0_gen = sparse(vcat(x_ind_A, x_ind_B), vcat(y_ind_A, y_ind_B), vcat(G0_values, G0_nulls))
            G1_gen = sparse(vcat(x_ind_A, x_ind_B), vcat(y_ind_B, y_ind_B), vcat(G1_values, G0_nulls))
            GAdd_gen = sparse(vcat(x_ind_A, x_ind_B), vcat(y_ind_A, y_ind_B), vcat(GAdd_values, G0_nulls))
            return G0_gen, G1_gen, GAdd_gen

        else
            G0_gen = G0
            G1_gen = G1
            GAdd_gen = GAdd

            return G0_gen, G1_gen, GAdd_gen
        end
    end
    G0_wth_div, G1_wth_div, GAdd_wth_div = G_gen(div)
    # dropzeros!(G0_wth_div;trim=true)
    # dropzeros!(G1_wth_div;trim=true)
    # dropzeros!(GAdd_wth_div;trim=true)
    return G0_wth_div, G1_wth_div, GAdd_wth_div
end