# This file contains code that was formerly a part of Julia.
# License is MIT: https://julialang.org/license

module KahanSummation

export sum_kbn, cumsum_kbn

"""
    cumsum_kbn(A; dims=:)

Cumulative sum, optionally along a dimension, using the Kahan-Babuska-Neumaier compensated
summation algorithm for additional accuracy.
"""
function cumsum_kbn end

cumsum_kbn(x::AbstractArray; dims=:) = _cumsum_kbn(x, dims)
cumsum_kbn(x; dims=:) = _cumsum_kbn(collect(x), dims)

function _cumsum_kbn(A::AbstractArray{T}, dims::Integer) where {T}
    axis_size = size(A, dims)
    axis_stride = 1
    for i = 1:dims-1
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

function _cumsum_kbn(v::AbstractArray{T}, ::Colon) where {T}
    r = similar(v)
    isempty(v) && return r
    inds = axes(v, 1)
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

"""
    TwicePrecisionN{T}(number)
    TwicePrecisionN{T}(hi, nlo)

Represents an extended precision number as `x.hi - x.nlo`.
We store the lower order component as the negation to avoid problems when `x.hi == -0.0`.

This does not subtype Number or Real, being meant primarily for internal use
in KahanSummation.jl.

Convert a `TwicePrecisionN{T}` back to a `T` by calling `singleprec(tp)`.
"""
struct TwicePrecisionN{T}
    hi::T
    nlo::T
end

singleprec(x::TwicePrecisionN{T}) where {T} = convert(T, x)

# Implement Base methods
Base.convert(::Type{TwicePrecisionN{T}}, x::Number) where {T} =
    TwicePrecisionN{T}(convert(T, x), zero(T))
Base.convert(::Type{T}, x::TwicePrecisionN) where {T} =
    convert(T, x.hi - x.nlo)

# Two-sum implementation
@inline function plus_kbn(x::T, y::T) where {T}
    hi = x + y
    nlo = abs(x) > abs(y) ? (hi - x ) - y : (hi - y) - x
    TwicePrecisionN(hi, nlo)
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

# Implement methods for accumulators, specifically mapreduce
Base.mapreduce_empty(f, ::typeof(plus_kbn), T) = TwicePrecisionN(zero(T),zero(T))
Base.mapreduce_empty(::typeof(identity), ::typeof(plus_kbn), T) = TwicePrecisionN(zero(T),zero(T)) # disambiguate
Base.mapreduce_first(f, ::typeof(plus_kbn), x) = TwicePrecisionN(x, zero(x))

# Finally, the implementation of `sum_kbn` is trivial, dispatching to `mapreduce`.  
# Most of the work happens in `plus_kbn`.

"""
    sum_kbn([f,] A)

Return the sum of all elements of `A`, using the Kahan-Babuska-Neumaier compensated
summation algorithm for additional accuracy.
"""
sum_kbn(f, X; kw..) = singleprec(mapreduce(f, plus_kbn, X; kw...))
sum_kbn(X; kw...) = sum_kbn(identity, X; kw...)


### Deprecations

Base.@deprecate cumsum_kbn(A, axis) cumsum_kbn(A; dims=axis)

end # module
