function additiveFromStoch(G, a, w)
    itoisometrymethod = Trapezoidal(innerSubdiv, 1.0)

    mSize = size(G)[1]
    M = zeros(mSize, mSize)
    for i = 1:mSize
        for j = 1:mSize
            try
                M[i, j] = sum([G[i, k].idx == w[j, j].idx ? itoisometrymethod(G[i, k], w[j, j]) * a[k] : 0 for k = 1:mSize])
            catch
                (err)

            end
        end
    end

    return M + M'
end

function secondMoment(coeffM::CoeffMatrices, collMethod, order, m=0; calc_time=false, kwargs...)
    # println(coeffM.G1)
    F0_map, F1_map, cMat_map = sparse.(M2_Mapping_from_Sparse.(sparse.([coeffM.F0, coeffM.F1, coeffM.cMat])))
    G0_map, G1_map, σMat_map = sparse.(M2_Mapping_from_Sparse.([coeffM.G0, coeffM.G1, coeffM.σMat]))
    nodes = collMethod.nodes
    divs = collMethod.div
    sizeM = nodes * divs * order
    e = zeros(sizeM, 1)
    if m == 0
        # we dont handle periodic solution
        e[1:nodes*order, 1] = ones(nodes * order, 1)
    else
        # we handle periodic solution, because "m" is given
        sizeM = sizeM * m + (divs - 1) * nodes * order
        e = ones(sizeM, 1)
    end
    σVec = σMat_map * mapMatToVec(e * e')
    cVec = cMat_map * mapMatToVec(e * e')

    H0 = Matrix(F0_map - G0_map)
    H1 = Matrix(F1_map + G1_map)
    H = H0 \ H1
    ρ = 0
    try
        ρ = abs(eigs(H, nev=1)[1][1])
    catch
        ρ = maximum(abs.(LinearAlgebra.eigvals(H)))
    end

    # time the eigenvalue calculation
    time = 0.0
    if calc_time
        time = @elapsed (abs(eigs(H, nev=1)[1][1]))
    end

    FirstM = firstMoment(coeffM, collMethod, order, m)
    FM = FirstM[2]
    solFM(x) = FirstM[3](x)
    additiveFromFM = mapMatToVec(coeffM.F1 * FM * diag(coeffM.cMat)' + diag(coeffM.cMat) * FM' * coeffM.F1')
    additive(x) = mapMatToVec(coeffM.F1 * solFM(x) * diag(coeffM.cMat)' + diag(coeffM.cMat) * solFM(x)' * coeffM.F1')

    # additiveStochMapping = mapMatToVec(additiveFromStoch(coeffM.G1, FM, coeffM.σMat)) + mapMatToVec(additiveFromStoch(coeffM.G0, FM, coeffM.σMat));
    additiveStochMapping = mapMatToVec(additiveFromStoch(coeffM.G1, FM, coeffM.σMat))

    additiveSt(x) = mapMatToVec(additiveFromStoch(coeffM.G1, solFM(x), coeffM.σMat))
    # G01 = M2_Mapping_from_Sparse(coeffM.G0, coeffM.G1)
    # H01 = G01 + G01'
    sol(x2, x) = H0 \ (H1 * mapMatToVec(x2) + cVec + σVec + additive(x) + additiveSt(x))
    # the (H01 + H01') * mapMatToVec(x2) part is only an approximation, bc a0*a1 is present

    # statSol = (H0 - H01 - H1) \ (cVec + σVec + additiveFromFM + additiveStochMapping);
    statSol = (H0 - H1) \ (cVec + σVec + additiveFromFM + additiveStochMapping)

    # return ρ, MapVecToMat(statSol)[1:order:end,1:order:end], MapVecToMat(statSol)[2:order:end,2,order:end]
    return ρ, statSol, sol, additiveSt, time
end