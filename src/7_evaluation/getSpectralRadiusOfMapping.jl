function getSpectralRadiusOfMapping(fm, sm)
    return (
        ρ_fm=fm[1],
        ρ_sm=sm[1]
    )
end

function getSpectralRadiusOfMapping(sm::SMMapping)
    # present
    H = Matrix(sm.F0 - sm.G0)
    # delay
    Hτ = Matrix(sm.F1 + sm.G1)
    # mapping matrix
    H2 = H \ Hτ
    ρ = 0.0
    try
        ρ = abs(eigs(H2, nev=1)[1][1])
    catch
        ρ = maximum(abs.(LinearAlgebra.eigvals(H2)))
    end
    return ρ
end

function getSpectralRadiusOfMapping(fm::FMMapping)
    H1 = fm.F0 \ fm.F1
    ρ = 0
    try
        ρ = abs(eigs(H1, nev=1)[1][1])
    catch
        ρ = maximum(abs.(LinearAlgebra.eigvals(H1)))
    end
    return ρ
end