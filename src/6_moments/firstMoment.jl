function firstMoment(coeffM::CoeffMatrices, collMethod, order, m=0)
    H = coeffM.F0 \ coeffM.F1
    ρ = 0
    try
        ρ = abs(eigs(H, nev=1)[1][1])
    catch
        ρ = maximum(abs.(LinearAlgebra.eigvals(H)))
    end
    # nodes = collMethod.nodes;
    # divs = collMethod.div;
    # sizeM = nodes*divs*order;
    # if m>0
    #     sizeM = sizeM*m+(divs-1)*nodes*order
    # end
    # e = ones(sizeM,1);
    # cVec = coeffM.cMat*e;
    cVec = diag(coeffM.cMat)
    statSol = (coeffM.F0 - coeffM.F1) \ (cVec)
    sol(x) = coeffM.F0 \ (coeffM.F1 * x + cVec)
    return ρ, statSol, sol
end