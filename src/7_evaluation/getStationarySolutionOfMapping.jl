struct STSolution{a_stT,a_fT,x_stT,x_stfT,a_st_HT,a_st_fT}
    a_st::a_stT
    a_f::a_fT
    x_st::x_stT
    x_stf::x_stfT
    a_st_H::a_st_HT
    a_st_f::a_st_fT
end

function additiveFromStoch(G::GT, a::aT, w::wT) where {GT,aT,wT}
    itoisometrymethod = Trapezoidal(innerSubdiv, 1.0)

    s = size(G)[1]
    M = zeros(s, s)
    for i = 1:s
        for j = 1:s
            try
                M[i, j] = sum([G[i, k].idx == w[j, j].idx ? itoisometrymethod(G[i, k], w[j, j]) * a[k] : 0 for k = 1:s])
            catch
                (err)

            end
        end
    end

    return mapMatToVec(M + M')
end

function additiveFromFM(F::FT, FM::FMT, C::CT) where {FT,FMT,CT}
    return mapMatToVec(F * FM * diag(C)' + diag(C) * FM' * F')
end

function getStationarySolutionOfMapping(smm::SMMapping, fmm::FMMapping, CM::CollocationMethod{xsT,lsT}, order, m=0; use_n=undef, kwargs...) where {xsT,lsT}
    s = CM.nodes * CM.div * order
    e = zeros(s, 1)
    if m == 0
        # no periodic solution
        e[1:CM.nodes*order, 1] = ones(CM.nodes * order, 1)
    else
        s = s * m + (CM.div - 1) * CM.nodes * order
        e = ones(sie, 1)
    end
    ev = mapMatToVec(e * e')
    w = smm.w * ev
    f = smm.f * ev

    fm_st = getStationarySolutionOfMapping(fmm, CM, order, m)
    add_fm = additiveFromFM(fmm.F1, fm_st.a_st, fmm.f)
    add_fm_f(x) = additiveFromFM(fmm.F1, fm_st.a_f(x), fmm.f)

    add_stoch = additiveFromStoch(fmm.G1, fm_st.a_st, fmm.w)
    add_stoch_f(x) = additiveFromStoch(fmm.G1, fm_st.a_f(x), fmm.w)

    H = Matrix(smm.F0 - smm.G0)
    Hτ = Matrix(smm.F1 + smm.G1)

    # nodal solutions
    # for testing----
    a_st_H = H - Hτ
    a_st_f = f + w + add_fm + add_stoch
    # -------
    a_st = zeros(size(f))
    if use_n == undef
        a_st = (H - Hτ) \ (f + w + add_fm + add_stoch)
    else
        a_st[1:use_n] = (H-Hτ)[1:use_n, 1:use_n] \ (f+w+add_fm+add_stoch)[1:use_n]
    end
    a_f(aτ2, aτ) = H \ (Hτ * mapMatToVec(aτ2) + f + w + add_fm_f(aτ) + add_stoch_f(aτ))

    ϕ(x) = [ϕi(x) for ϕi in CM.ls]

    # approximated state solutions
    a_st_m = mapVecToMat(a_st)
    x_stf(x1, x2) = [ϕ(x1)' * a_st_m[o:order:end, o:order:end] * ϕ(x2) for o in 1:order]
    x_st = [x_stf(x1, x2) for x1 in CM.xs, x2 in CM.xs]

    return STSolution(
        a_st, a_f, x_st, x_stf, a_st_H, a_st_f
    )

end

function getStationarySolutionOfMapping(fmm::FMMapping, CM::CollocationMethod{xsT,lsT}, order, m=0) where {xsT,lsT}
    f = diag(fmm.f)
    # nodal solutions
    a_st_H = fmm.F0 - fmm.F1
    a_st_f = f
    a_st = (fmm.F0 - fmm.F1) \ f
    a_f(aτ) = fmm.F0 \ (fmm.F1 * aτ + f)

    ϕ(x) = [ϕi(x) for ϕi in CM.ls]

    # state approximation solutions
    x_stf(x) = [sum(ϕ(x)'a_st[o:order:end]) for o in 1:order]
    x_st = x_stf.(CM.xs)

    return STSolution(
        a_st, a_f, x_st, x_stf, a_st_H, a_st_f
    )
end