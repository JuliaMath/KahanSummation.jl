# This file contains code that was formerly a part of Julia.
# License is MIT: https://julialang.org/license

__precompile__(true)

module KahanSummation

if VERSION >= v"0.7.0-DEV.3000" # TODO: More specific bound
    if isdefined(Base, :sum_kbn) # Deprecated
        import Base: sum_kbn, cumsum_kbn
    else
        export sum_kbn, cumsum_kbn
    end
end

if isdefined(Base, Symbol("@default_eltype"))
    using Base: @default_eltype
else
    macro default_eltype(itr)
        quote
            Core.Inference.return_type(first, Tuple{$(esc(itr))})
        end
    end
end

if isdefined(Base, :promote_sys_size_add)
    using Base: promote_sys_size_add
else
    promote_sys_size_add(x::T) where {T} = Base.r_promote(+, zero(T)::T)
end

"""
    TwicePrecisionN{T}

Represents an extended precision number as `x.hi - x.nlo`.

We store the lower order component as the negation to avoid problems when `x.hi == -0.0`.
"""
struct TwicePrecisionN{T}
    hi::T
    nlo::T
end


@inline function plus_kbn(x::T, y::T) where {T}
    hi = x + y
    nlo = abs(x) > abs(y) ? (hi - x ) - y : (hi - y) - x
    TwicePrecisionN(hi, lo)
end
@inline function plus_kbn(x::T, y::TwicePrecisionN{T}) where {T}
    hi = x + y.hi
    if abs(x) > abs(y.hi)
        nlo = ((hi - x) - y.hi) + y.nlo
    else
        nlo = ((hi - y.hi) - x) + y.nlo
    end
    TwicePrecisionN(hi, nlo)
end
@inline plus_kbn(x::TwicePrecisionN{T}, y::T) where {T} = plus_kbn(y, x)

@inline function plus_kbn(x::TwicePrecisionN{T}, y::TwicePrecisionN{T}) where {T}
    hi = x.hi + y.hi
    if abs(x.hi) > abs(y.hi)
        nlo = (((hi - x.hi) - y.hi) + y.nlo) + x.nlo
    else
        nlo = (((hi - y.hi) - x.hi) + x.nlo) + y.nlo
    end
    TwicePrecisionN(hi, nlo)
end

Base.r_promote_type(::typeof(plus_kbn), ::Type{T}) where {T} =
    TwicePrecisionN{T}

Base.convert(::Type{TwicePrecisionN{T}}, x::Number) where {T} =
    TwicePrecisionN{T}(convert(T, x), zero(T))
Base.convert(::Type{T}, x::TwicePrecisionN) where {T} =
    convert(T, x.hi - x.nlo)

Base.mr_empty(f, ::typeof(plus_kbn), T) = TwicePrecisionN(zero(T),zero(T))

singleprec(x::TwicePrecisionN{T}) where {T} = convert(T, x)


"""
    sum_kbn([f,] A)

Return the sum of all elements of `A`, using the Kahan-Babuska-Neumaier compensated
summation algorithm for additional accuracy.
"""
sum_kbn(f, X) = singleprec(mapreduce(f, plus_kbn, X))
sum_kbn(X) = sum_kbn(identity, X)







"""
    cumsum_kbn(A, dim::Integer)

Cumulative sum along a dimension, using the Kahan-Babuska-Neumaier compensated summation
algorithm for additional accuracy.
"""
function cumsum_kbn(A::AbstractArray{T}, axis::Integer) where T<:AbstractFloat
    dimsA = size(A)
    ndimsA = ndims(A)
    axis_size = dimsA[axis]
    axis_stride = 1
    for i = 1:(axis-1)
        axis_stride *= size(A, i)
    end
    axis_size <= 1 && return A
    B = similar(A)
    C = similar(A)
    for i = 1:length(A)
        if div(i-1, axis_stride) % axis_size == 0
            B[i] = A[i]
            C[i] = zero(T)
        else
            s = B[i-axis_stride]
            Ai = A[i]
            B[i] = t = s + Ai
            if abs(s) >= abs(Ai)
                C[i] = C[i-axis_stride] + ((s-t) + Ai)
            else
                C[i] = C[i-axis_stride] + ((Ai-t) + s)
            end
        end
    end
    return B + C
end

function cumsum_kbn(v::AbstractVector{T}) where T<:AbstractFloat
    r = similar(v)
    isempty(v) && return r
    inds = indices(v, 1)
    i1 = first(inds)
    s = r[i1] = v[i1]
    c = zero(T)
    for i = i1+1:last(inds)
        vi = v[i]
        t = s + vi
        if abs(s) >= abs(vi)
            c += ((s-t) + vi)
        else
            c += ((vi-t) + s)
        end
        s = t
        r[i] = s+c
    end
    return r
end

end # module
