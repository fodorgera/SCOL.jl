"""
`stabilityMapFunction(paramX, paramY ,system, params, coll)`

* `paramX` and `paramY` are used as the parameters of the stability map
* `params` contain all the parameters, not only the remaining ones
"""
function stabilityMapFunction(paramX, paramY, system, params, coll)
    # [2024-02-27] - function retructures -> does not work anymore
    mapping = getMappingOfSystem(system, [paramX, paramY, params[3:end]...], coll)
    spectralRadius = getSpectralRadiusOfMapping(mapping...)
    ρ_sm = spectralRadius.ρ_sm
    return ρ_sm
end

function getStabilityMapOfSystem(foo::Function; axis = [Axis(-15.001:15.001, :paramX), Axis(-15.001:15.001, :paramY)], iteration = 4, kwargs...)
    x, y = getinterpolatedsolution(solve!(MDBM_Problem(foo, axis), iteration))
    return x, y
end