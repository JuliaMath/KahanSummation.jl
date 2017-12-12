# This file contains code that was formerly a part of Julia.
# License is MIT: https://julialang.org/license

using KahanSummation
using Compat
using Compat.Test

# Shadow the names in Base if they're defined
import KahanSummation: sum_kbn, cumsum_kbn

@testset "cumsum_kbn" begin
    v   = [1,1e100,1,-1e100]*1000
    v2  = [1,-1e100,1,1e100]*1000

    cv  = [1,1e100,1e100,2]*1000
    cv2 = [1,-1e100,-1e100,2]*1000

    @test isequal(cumsum_kbn(v), cv)
    @test isequal(cumsum_kbn(v2), cv2)

    A = [v reverse(v) v2 reverse(v2)]

    c = cumsum_kbn(A, 1)

    @test isequal(c[:,1], cv)
    @test isequal(c[:,3], cv2)
    @test isequal(c[4,:], [2.0, 2.0, 2.0, 2.0]*1000)

    c = cumsum_kbn(A, 2)

    @test isequal(c[1,:], cv2)
    @test isequal(c[3,:], cv)
    @test isequal(c[:,4], [2.0,2.0,2.0,2.0]*1000)
    
    amat = reshape(collect(1.0:16.0),(4,4))
  
    @test cumsum_kbn(amat, 2)[3,:] == [3.0, 10.0, 21.0, 36.0]
    @test diag(cumsum_kbn(amat, 1)) == [1.0, 11.0, 30.0, 58.0]

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

    @test sum_kbn(identity, [1,1e100,1,-1e100]) === 2.0
    @test sum_kbn(identity, Float64[]) === 0.0
    @test sum_kbn(identity, i for i=1.0:1.0:10.0) === 55.0
    @test sum_kbn(identity, i for i=1:1:10) === 55
    @test sum_kbn(identity, [1 2 3]) === 6
    @test sum_kbn(identity, [2+im 3-im]) === 5+0im
    @test sum_kbn(identity, [1+im 2+3im]) === 3+4im
    @test sum_kbn(identity, [7 8 9]) === sum_kbn([8 9 7])
    @test sum_kbn(identity, i for i=1:1:10) === sum_kbn(i for i=10:-1:1)
    @test sum_kbn(identity, [-0.0]) === -0.0
   
    @test sum_kbn([-0.0]) === -0.0
    @test sum_kbn(x->x*x, [1.0]) === 1.0
    @test sum_kbn(x->x*x, [2.0]) === 4.0
    @test sum_kbn(x->x+1, Float64[]) === 0.0
    @test sum_kbn(x->x+1, [6 7 8]) === sum_kbn([8 9 7])
    @test sum_kbn(sqrt, i for i=1.0:1.0:10.0) === 22.4682781862041
    @test sum_kbn(x->x*im,[1+im 2+3im]) === -4+3im
    @test sum_kbn(sqrt, [i for i in 1.0:1.0:10.0]) === sum_kbn(sqrt, [i for i in 10.0:-1.0:1.0])
    @test sum_kbn([1.0, 1.0e16, 1.0, -1.0e16]) === sum_kbn(x->x+1, [1.0, 1.0e16, 1.0, -1.0e16] - 1.0)
end
