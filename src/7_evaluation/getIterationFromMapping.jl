function getIterationFromMapping(fm, sm, Coll::CollocationMethod; segments=5, a0=[idx == 1 ? 1.0 : 0.0 for idx in eachindex(fm[2])], order=2, kwargs...)
    nodal_sols = [a0 for i in 1:segments]
    nodal_cov_sols = [a0 * a0' for i in 1:segments]

    for i in 2:segments
        nodal_sols[i] = fm[3](nodal_sols[i-1])
        nodal_cov_sols[i] = mapVecToMat(sm[3](nodal_cov_sols[i-1], nodal_sols[i-1]))
        nodal_cov_sols[i][1:order, 1:order] += nodal_cov_sols[i-1][end-order+1:end, end-order+1:end]
    end

    ϕs(x) = [li(x) for li in Coll.ls]
    sim_fm_st = x -> [sum(fm[2][i:order:end] .* ϕs(x)) for i in 1:order]
    sim_sm_st = x -> [ϕs(x)' * mapVecToMat(sm[2])[i:order:end, i:order:end] * ϕs(x) for i in 1:order]

    sim_fm_1 = [x -> sum(nodal_sols[i][1:order:end] .* ϕs(x)) for i in 1:segments]
    sim_sm_1 = [x -> ϕs(x)' * nodal_cov_sols[i][1:order:end, 1:order:end] * ϕs(x) for i in 1:segments]
    stD_1 = [x -> sqrt(abs(sim_sm_1[i](x)[1:order:end] - sim_fm_1[i](x)[1:order:end]^2)) for i in 1:segments]
    stD_st = x -> sqrt(abs(sim_sm_st(x)[1:order:end] - sim_fm_st(x)[1:order:end]^2))

    return sim_fm_1, sim_sm_1, stD_1, sim_fm_st, sim_sm_st, stD_st
end