function buildStochMXsWholeState(αt::Function, βt::Function, σt::Function, m, τ, Coll)::StochMXs
    # ToDo x,y has to be shifted according to j: (j-1)*nodes*order*div
    αtMX = vcat([buildSparseFromFunction2(αt, j, τ, Coll) for j = 1:m]...)
    βtMX = vcat([buildSparseFromFunction2(βt, j, τ, Coll) for j = 1:m]...)
    σtVEC = vcat([buildAddVecFromFunction(σt, j, τ, Coll) for j = 1:m]...)
    # Constructing the indices
    order = size(αt(0.0))[1]
    x_ind, y_ind = sparseIndexes(order, Coll.nodes)
    x = vcat([x_ind .+ (Coll.nodes * order * (k - 1)) for k = 1:Coll.div]...)
    y = vcat([y_ind .+ (Coll.nodes * order * (k - 1)) for k = 1:Coll.div]...)

    # updating the x,y indices according to j
    xs = vcat([x .+ ((j - 1) * Coll.nodes * order * Coll.div * order) for j = 1:m]...)
    ys = vcat([y .+ ((j - 1) * Coll.nodes * order * Coll.div * order) for j = 1:m]...)
    # shifting according to the whole state
    xsα = xs .+ (Coll.div - 1) * Coll.nodes * order
    ysα = ys .+ (Coll.div - 1) * Coll.nodes * order
    xsβ = xs .+ (Coll.div - 1) * Coll.nodes * order
    ysβ = ys
    if Coll.div > 1
        append!(xsβ, Coll.nodes * (2 * Coll.div - 1) * order + (m - 1) * Coll.nodes * Coll.div * order)
        append!(ysβ, Coll.nodes * (2 * Coll.div - 1) * order + (m - 1) * Coll.nodes * Coll.div * order)
        subdiv = innerSubdiv
        βMXadd = ItoSVec(0, SVector{subdiv + 1}([0.0 for i = 1:subdiv+1]))
        return StochMXs(sparse(xsα, ysα, αtMX), sparse(xsβ, ysβ, [βtMX..., βMXadd]), sparse(xsα, ysα, σtVEC))
    else
        return StochMXs(sparse(xsα, ysα, αtMX), sparse(xsβ, ysβ, βtMX), sparse(xsα, ysα, σtVEC))
    end
end