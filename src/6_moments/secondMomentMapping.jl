struct SMMapping{F0T,F1T,fT,G0T,G1T,wT}
    F0::F0T
    F1::F1T
    f::fT
    G0::G0T
    G1::G1T
    w::wT
end

function secondMomentMapping(coeffM::CoeffMatrices)
    F0_map, F1_map, cMat_map = dropzeros.(M2_Mapping_from_Sparse.(sparse.([coeffM.F0, coeffM.F1, coeffM.cMat])))
    G0_map, G1_map, σMat_map = dropzeros.(M2_Mapping_from_Sparse.([coeffM.G0, coeffM.G1, coeffM.σMat]))
    return SMMapping(
        F0_map, F1_map, cMat_map,
        G0_map, G1_map, σMat_map
    )
end