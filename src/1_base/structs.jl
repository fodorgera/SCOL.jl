struct CollocationMethod{xsT,lsT}
    name::String
    nodes::Int64
    div::Int64

    xs::xsT
    ls::lsT
    mult::Float64

    MXPack::Tuple{Array{Float64,2},Array{Float64,2},Array{Float64,2}}
end

struct SCOLProblem{AT,BT,cT,αT,βT,γT,τT,TT}
    A::AT
    B::BT
    c::cT
    α::αT
    β::βT
    γ::γT
    τ::τT
    T::TT
end