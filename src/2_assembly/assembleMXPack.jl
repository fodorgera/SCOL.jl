function assembleMXPack(nodes, xs, ls) where {xsT,lsT}
    MX_func = zeros(nodes, nodes)
    MX_end = copy(MX_func)
    MX_start = copy(MX_func)
    for i = 1:nodes-1
        MX_func[i+1, :] = [quadgk(x -> ls[j](x), xs[i], xs[i+1], atol = 0.0001)[1] for j = 1:nodes]
        MX_end[i+1, :] = [ls[j](xs[i+1]) for j = 1:nodes]
        MX_start[i+1, :] = [ls[j](xs[i]) for j = 1:nodes]
    end
    MX_end[1, :] += xs[1] .|> ls
    return MX_end, MX_start, MX_func
end