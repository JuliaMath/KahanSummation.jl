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
    sum_kbn(A)

Return the sum of all elements of `A`, using the Kahan-Babuska-Neumaier compensated
summation algorithm for additional accuracy.
"""
function sum_kbn(A)
    T = Base.@default_eltype(A)
    c = Base.reduce_empty(+, T)
    it = iterate(A)
    it === nothing && return c
    Ai, i = it
    s = Ai - c
    while (it = iterate(A, i)) !== nothing
        Ai, i = it::Tuple{T, Int}
        t = s + Ai
        if abs(s) >= abs(Ai)
            c -= ((s-t) + Ai)
        else
            c -= ((Ai-t) + s)
        end
        s = t
    end
    s - c
end

### Deprecations

Base.@deprecate cumsum_kbn(A, axis) cumsum_kbn(A; dims=axis)

end # module
