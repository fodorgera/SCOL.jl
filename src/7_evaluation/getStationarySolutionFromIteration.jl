function getStationarySolutionFromIteration(iteration)
    fm_st = iteration[4]
    sm_st = iteration[5]
    stD_st = iteration[6]
    return (
        fm_st = fm_st,
        sm_st = sm_st,
        stD_st = stD_st
    )
end