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

"""
    sum_kbn([f,] A)

Return the sum of all elements of `A`, using the Kahan-Babuska-Neumaier compensated
summation algorithm for additional accuracy.  When a function `f` is supplied, the
result of calling f on each element of `A` is summed.
"""
function sum_kbn(A)
    T = @default_eltype(typeof(A))
    c = promote_sys_size_add(zero(T)::T)
    i = start(A)
    if done(A, i)
        return c
    end
    Ai, i = next(A, i)
    s = Ai - c
    while !(done(A, i))
        Ai, i = next(A, i)
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

function sum_kbn(f,A)
    T = @default_eltype(typeof(A))
    c = promote_sys_size_add(zero(T)::T)
    i = start(A)
    if done(A, i)
        return c
    end
    Ai, i = next(A, i)
    s = f(Ai) - c
    while !(done(A, i))
        Ai, i = next(A, i)
        fAi = f(Ai); t = s + fAi
        if abs(s) >= abs(fAi)
            c -= ((s-t) + fAi)
        else
            c -= ((fAi-t) + s)
        end
        s = t
    end
    s - c
end

sum_kbn(identity, A) = sum_kbn(A)

end # module
