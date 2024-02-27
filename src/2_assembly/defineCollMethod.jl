function defineCollMethod(node_num, div_num, method_name)
    nodes = node_num
    div = div_num
    name = method_name

    mult, method = selectMethod(name, nodes)
    xs, ls = method
    MXPack = assembleMXPack(nodes, xs, ls)
    collMethod = CollocationMethod{typeof(xs),typeof(ls)}(name, nodes, div, xs, ls, mult, MXPack)
    return collMethod
end