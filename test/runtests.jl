# This file contains code that was formerly a part of Julia.
# License is MIT: https://julialang.org/license

using KahanSummation
using Test

@testset "cumsum_kbn" begin
    v   = [1,1e100,1,-1e100]*1000
    v2  = [1,-1e100,1,1e100]*1000

    cv  = [1,1e100,1e100,2]*1000
    cv2 = [1,-1e100,-1e100,2]*1000

    @test isequal(cumsum_kbn(v), cv)
    @test isequal(cumsum_kbn(v2), cv2)

    A = [v reverse(v) v2 reverse(v2)]

    c = cumsum_kbn(A; dims=1)

    @test isequal(c[:,1], cv)
    @test isequal(c[:,3], cv2)
    @test isequal(c[4,:], [2.0, 2.0, 2.0, 2.0]*1000)

    c = cumsum_kbn(A; dims=2)

    @test isequal(c[1,:], cv2)
    @test isequal(c[3,:], cv)
    @test isequal(c[:,4], [2.0,2.0,2.0,2.0]*1000)
    
    @test isequal(cumsum_kbn(1:3), [1,3,6])
    @test isequal(cumsum_kbn((i for i in [1,2,3])), [1,3,6])
    
end

@testset "sum_kbn" begin
    @test sum_kbn([1,1e100,1,-1e100]) === 2.0
    @test sum_kbn(Float64[]) === 0.0
    @test sum_kbn(i for i=1.0:1.0:10.0) === 55.0
    @test sum_kbn(i for i=1:1:10) === 55
    @test sum_kbn([1 2 3]) === 6
    @test sum_kbn([2+im 3-im]) === 5+0im
    @test sum_kbn([1+im 2+3im]) === 3+4im
    @test sum_kbn([7 8 9]) === sum_kbn([8 9 7])
    @test sum_kbn(i for i=1:1:10) === sum_kbn(i for i=10:-1:1)
    @test sum_kbn([-0.0]) === -0.0
    @test sum_kbn([-0.0,-0.0]) === -0.0
    @test sum_kbn(Iterators.filter(isodd, 1:10)) == 25
    @test isequal(sum_kbn(1:3), 6)
    @test isequal(sum_kbn((i for i in [1,2,3])), 6)
end

@testset "twice-precision addition"
    # Note: the functions and types used here are internal

    # The intent of this test is to make sure that the two-sum
    # works on two twiceprecision numbers correctly, nothing else.
    i1 = convert(KahanSummation.TwicePrecisionN{Float64}, 1e100)
    i2 = KahanSummation.TwicePrecisionN{Float64}(-1e100, 1)
    i12 = KahanSummation.plus_kbn(i1, i2)
    f12 = KahanSummation.singleprec(i12)
    @test f12 == 1
end