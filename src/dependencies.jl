using SemiDiscretizationMethod

innerSubdiv = 10

mutable struct ItoSVec{T,TV} <: AbstractVector{TV}
    idx::UInt64
    v::T
end

import Base: size, getindex
ItoSVec(idx, v) = ItoSVec{typeof(v),eltype(v)}(idx, v)
size(v::ItoSVec) = size(v.v)
getindex(v::ItoSVec, idx) = v.v[idx]

itoisometry_0(a::ItoSVec{T,TV}, b::ItoSVec{T,TV}) where {T,TV} = a.idx == b.idx ? Trapezoidal(innerSubdiv, 1.0)(a.v, b.v) : zero(TV)
#########################################

#########################################
# Funtions from M2_Mapping.jl
#########################################
abstract type ItoIsometryMethod{K} <: AbstractVector{Float64} end
struct Trapezoidal{K} <: ItoIsometryMethod{K}
    dt::Float64
end
Trapezoidal(k, T) = Trapezoidal{k + 1}(T / k)
Base.size(iim::Trapezoidal{k}) where {k} = (k,)
Base.getindex(iim::Trapezoidal{k}, idx::Integer) where {k} = idx < k && idx > 1 ? 1.0 : 0.5
Base.getindex(iim::Trapezoidal, idxs::Vector{<:Integer}) = getindex.(Ref(iim), idxs)
for T in union(InteractiveUtils.subtypes(ItoIsometryMethod))
    @eval (::$T)(val1::Real, val2::Real) = val1 * val2
end
# function (iim::Trapezoidal)(f1::AbstractVector{<:Real},f2::AbstractVector{<:Real})
#     sum(f1 .* f2 .* iim) * iim.dt
# end
function (iim::Trapezoidal)(f1::AbstractVector{<:Real}, f2::AbstractVector{<:Real})
    sum(f1 .* f2 .* iim) * iim.dt
end

struct CovVecIdx
    sectionStarts::Vector{Int64} # Every diagonal's start-1 in the covariance vector
    function CovVecIdx(siz::Int64)
        resvec = Vector{Int64}(undef, siz + 1)
        resvec[1] = 0
        for (i, n) in enumerate(siz:-1:1)
            resvec[i+1] = resvec[i] + n
        end
        new(resvec)
    end
    CovVecIdx(MX::AbstractArray) = CovVecIdx(size(MX, 1))
end

function (cvIdx::CovVecIdx)(i::Int64, j::Int64)
    cvIdx.sectionStarts[abs(i - j)+1] + min(i, j)
end

function M2_Mapping_from_Sparse(SM::SparseMatrixCSC{TSMX,Int64}, itoisometrymethod=Trapezoidal(innerSubdiv, 1.0)) where {TSMX<:Union{<:AbstractVector{T},<:T}} where {T<:Real}
    idx = CovVecIdx(SM)
    (_Is, _Js, _Vs) = findnz(SM)
    elemnum = sum(eachindex(_Vs))
    Is = Vector{Int64}(undef, elemnum)
    Js = Vector{Int64}(undef, elemnum)
    Vs = Vector{T}(undef, elemnum)
    k = 1
    for i in eachindex(_Is)
        for j = i:length(_Is)
            Is[k] = idx(_Is[i], _Is[j])
            Js[k] = idx(_Js[i], _Js[j])
            if Is[k] <= idx.sectionStarts[2] && Js[k] > idx.sectionStarts[2]
                Vs[k] = 2 * itoisometrymethod(_Vs[i], _Vs[j])
            else
                Vs[k] = itoisometrymethod(_Vs[i], _Vs[j])
            end
            k += 1
        end
    end
    sparse(Is, Js, Vs, idx.sectionStarts[end], idx.sectionStarts[end])
end

function M2_Mapping_from_Sparse(SM::SparseMatrixCSC{TSMX,Int64}, itoisometrymethod=itoisometry_0) where {TSMX<:Union{<:ItoSVec{T,TV}}} where {T,TV}
    idx = CovVecIdx(SM)
    (_Is, _Js, _Vs) = findnz(SM)
    elemnum = sum(eachindex(_Vs))
    Is = Vector{Int64}(undef, elemnum)
    Js = Vector{Int64}(undef, elemnum)
    Vs = Vector{TV}(undef, elemnum)
    k = 1
    for i in eachindex(_Is)
        for j = i:length(_Is)
            Is[k] = idx(_Is[i], _Is[j])
            Js[k] = idx(_Js[i], _Js[j])
            if Is[k] <= idx.sectionStarts[2] && Js[k] > idx.sectionStarts[2]
                Vs[k] = 2 * itoisometrymethod(_Vs[i], _Vs[j])
            else
                Vs[k] = itoisometrymethod(_Vs[i], _Vs[j])
            end
            k += 1
        end
    end
    sparse(Is, Js, Vs, idx.sectionStarts[end], idx.sectionStarts[end])
end

function M2_Mapping_from_Sparse(SM1::SparseMatrixCSC{TSMX,Int64}, SM2::SparseMatrixCSC{TSMX,Int64}, itoisometrymethod=itoisometry_0) where {TSMX<:Union{<:ItoSVec{T,TV}}} where {T,TV}
    idx = CovVecIdx(SM1)
    (_Is1, _Js1, _Vs1) = findnz(SM1)
    (_Is2, _Js2, _Vs2) = findnz(SM2)
    elemnum = sum(eachindex(_Vs1))
    Is = Vector{Int64}(undef, elemnum)
    Js = Vector{Int64}(undef, elemnum)
    Vs = Vector{TV}(undef, elemnum)
    k = 1
    for i in eachindex(_Is1)
        for j = i:length(_Is2)
            Is[k] = idx(_Is1[i], _Is1[j])
            Js[k] = idx(_Js1[i], _Js1[j])
            if Is[k] <= idx.sectionStarts[2] && Js[k] > idx.sectionStarts[2]
                Vs[k] = 2 * itoisometrymethod(_Vs1[i], _Vs2[j])
            else
                Vs[k] = itoisometrymethod(_Vs1[i], _Vs2[j])
            end
            k += 1
        end
    end
    sparse(Is, Js, Vs, idx.sectionStarts[end], idx.sectionStarts[end])
end