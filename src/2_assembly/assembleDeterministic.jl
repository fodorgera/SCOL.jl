function assembleDeterministic(A::Array{Float64,2}, B::Array{Float64,2}, c::Array{Float64}, τ::Float64, method::CollocationMethod{xsT,lsT}) where {xsT,lsT} # compatible for Chebyshev (hopefully)
    xs = method.xs
    ls = method.ls
    mult = method.mult
    nodes = method.nodes
    div = method.div

    # MX_end, MX_start, MX_func = MXassembly_func(method)
    MX_end = method.MXPack[1]
    MX_start = method.MXPack[2]
    MX_func = method.MXPack[3]
    diff_order = length(A[1, :]) # order of the differential equation

    Idiff = I + zeros(diff_order, diff_order)
    Idiv = I + zeros(div, div)

    # END POINTS MATRIX
    Φn = I + zeros(nodes * div, nodes * div)
    Φn[1:nodes, 1:nodes] = MX_end
    Φn_final = mxmult(Idiff, Φn)
    # START POINT MATRIX
    Φ0 = zeros(nodes * div, nodes * div)
    Φ0[1:nodes, 1:nodes] = MX_start
    Φ0_final = mxmult(Idiff, Φ0)
    # INTEGRATING THE CURRENT STATE
    AMX = zeros(nodes * div, nodes * div)
    AMX[1:nodes, 1:nodes] = τ * mult / div * MX_func
    AMX_final = mxmult(A, AMX)
    # INTEGRATING THE DELAY STATE
    BMX = zeros(nodes * div, nodes * div)
    BMX[1:nodes, end-nodes+1:end] = τ * mult / div * MX_func
    BMX_final = mxmult(B, BMX)
    # CONNECTION MATRIX
    CMX = zeros(nodes * div, nodes * div)
    CMX[1, 1:nodes] = xs[end] .|> ls
    CMX[nodes+1:end, 1:end-nodes] = I + zeros((div - 1) * nodes, (div - 1) * nodes)
    CMX_final = mxmult(Idiff, CMX)

    # ADDITIVE VECTOR
    cADD = zeros(nodes * div)
    cADD[2:nodes] = (τ * mult / div) * [xs[i] - xs[i-1] for i = 2:nodes]
    cADD_Mat = mxmult(diagm(c), diagm(cADD))
    # println(cADD_Mat)
    F0 = Φn_final - Φ0_final - AMX_final
    F1 = CMX_final + BMX_final

    return F0, F1, cADD_Mat
end