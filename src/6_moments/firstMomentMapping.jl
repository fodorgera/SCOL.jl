struct FMMapping{F0T,F1T,fT,G0T,G1T,wT}
    F0::F0T
    F1::F1T
    f::fT
    G0::G0T
    G1::G1T
    w::wT
end

function firstMomentMapping(coeffM::CoeffMatrices)
    return FMMapping(
        coeffM.F0, coeffM.F1, coeffM.cMat,
        coeffM.G0, coeffM.G1, coeffM.σMat
    )
end