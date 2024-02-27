function buildDetMXsWholeState(At::Function, Bt::Function, ct::Function, m, τ, Coll)::DetMXs
    nodes = Coll.nodes
    div = Coll.div
    order = size(At(0.0))[1]
    # considering the delay states as well in the state
    dim = nodes * div * order

    MX_end = Coll.MXPack[1]
    MX_start = Coll.MXPack[2]


    AtMX = [buildMXfromFunction(At, j, τ, Coll) for j = 1:m]
    BtMX = [buildMXfromFunction(Bt, j, τ, Coll) for j = 1:m]
    ctVEC = [buildVECfromFunction(ct, j, τ, Coll) for j = 1:m]

    # create the endpoint matrix
    EndMX = zeros(dim, dim)
    EndMX_k = mxmult(I + zeros(order, order), Coll.MXPack[1])
    ConnectionMX_k = zeros(nodes * order, nodes * order)
    ConnectionMX_k[1:order, :] = EndMX_k[end-order+1:end, :]
    # the whole detail about the connections
    ConnectionMX = mxmult(ConnectionMX_k, I + zeros(div, div))

    StartMX_k = mxmult(I + zeros(order, order), MX_start)
    StartMX_j = mxmult(StartMX_k, I + zeros(div, div))

    F0 = zeros(dim * m, dim * m)
    F1 = zeros(dim * m, dim * m)
    CVec = zeros(dim * m) # TODO: the 0s at the top has to be implemented
    for j = 1:m
        for k = 1:div
            # filling up the endpoint matrix
            EndMX[(k-1)*nodes*order+1:k*nodes*order, (k-1)*nodes*order+1:k*nodes*order] = EndMX_k
        end
        AtMX[j] += StartMX_j
        # filling up the large matrices
        F0[(j-1)*dim+1:j*dim, (j-1)*dim+1:j*dim] = EndMX - AtMX[j]
        F1[(j-1)*dim+1:j*dim, (j-1)*dim+1:j*dim] = BtMX[j]
        CVec[(j-1)*dim+1:j*dim] = ctVEC[j]
    end
    F0_full = I + zeros(nodes * (2 * div - 1) * order, nodes * (2 * div - 1) * order)
    F0_full[(div-1)*nodes*order+1:end, (div-1)*nodes*order+1:end] = F0
    F1_full = zeros(nodes * (2 * div - 1) * order, nodes * (2 * div - 1) * order)
    F1_full[(div-1)*nodes*order+1:end, (div-1)*nodes*order+1:end] += ConnectionMX
    F1_full[(div-1)*nodes*order+1:end, 1:end-(div-1)*nodes*order] += BtMX[1]
    F1_full[1:(div-1)*nodes*order, nodes*order+1:div*nodes*order] += I + zeros((div - 1) * nodes * order, (div - 1) * nodes * order)
    CVec_full = zeros(nodes * (2 * div - 1) * order)
    CVec_full[(div-1)*nodes*order+1:end] += CVec
    return DetMXs(F0_full, F1_full, diagm(CVec_full))
end