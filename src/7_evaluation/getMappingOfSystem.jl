struct SCOLMapping{fmT,smT}
    fm::fmT
    sm::smT
end

function getMappingOfSystem(sys::SCOLProblem,Coll::CollocationMethod{xsT,lsT}) where {xsT,lsT}
    try
        return getConstantMappingOfSystem(sys,Coll)
    catch err
        println("[SCOL][ERROR][getMappingOfSystem] using periodic (time dependent) solver.")
        return getPeriodicMappingOfSystem(sys,Coll)
    end
end

function getConstantMappingOfSystem(sys::SCOLProblem, Coll::CollocationMethod{xsT,lsT}; kwargs...) where {xsT,lsT}
    coeffM = buildCoeffMatrices(sys.A, sys.B, sys.c, sys.α, sys.β, sys.γ, sys.τ, Coll)
    # FirstM_P = firstMoment(coeffM, Coll, order, 1)
    mapping_fm = firstMomentMapping(coeffM)
    # SecondM_P = secondMoment(coeffM, Coll, order, 1, calc_time=calc_time)
    mapping_sm = secondMomentMapping(coeffM)
    return SCOLMapping(mapping_fm, mapping_sm)
    # return (
    #     mapping_fm=mapping_fm,
    #     mapping_sm=mapping_sm
    # )
end

function getPeriodicMappingOfSystem(sys::SCOLProblem, Coll::CollocationMethod{xsT,lsT}; kwargs...) where {xsT,lsT}
    # THERE IS SOMETHING WRONG WITH THE PERIODIC SOLVER -> SPARSE INDEXES, OR MX MULTIPLICATION?
    # 2022.07.09 - maybe it is okay no -> the x,y was reversed in the mxListIntoSparse

    DetMX = buildDetMXsWholeState(
        sys.A,
        sys.B,
        sys.c,
        1, sys.τ, Coll
    )
    StochMX = buildStochMXsWholeState(
        sys.α,
        sys.β,
        sys.γ,
        1, sys.τ, Coll
    )
    coeffM_P = CoeffMatrices(DetMX.F0, DetMX.F1, DetMX.ct, StochMX.αt, StochMX.βt, StochMX.σt)

    # coeffM = buildCoeffMatrices(sys.A, sys.B, sys.c, sys.α, sys.β, sys.γ, sys.τ, Coll)
    # FirstM_P = firstMoment(coeffM, Coll, order, 1)
    mapping_fm = firstMomentMapping(coeffM_P)
    # SecondM_P = secondMoment(coeffM, Coll, order, 1, calc_time=calc_time)
    mapping_sm = secondMomentMapping(coeffM_P)
    return SCOLMapping(mapping_fm, mapping_sm)
    # return (
    #     mapping_fm=mapping_fm,
    #     mapping_sm=mapping_sm
    # )
end