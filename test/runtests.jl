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
