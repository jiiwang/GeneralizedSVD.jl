include("../src/wrapper.jl")

using BenchmarkTools
using LinearAlgebra

function test(n)
    x = randn(n,1)[:,1];
    y = randn(n,1)[:,1];
    println(BLAS.dot(n,x,1,y,1) == compute_dot_julia(x, y))
    @btime compute_dot($x, $y)
    @btime BLAS.dot($n, $x, 1, $y, 1)
    @btime compute_dot_julia($x, $y)
end

function for_dot(x, y)
    dot = sum(x.*y)
end

function test1(n)
    maxJ = 20
    time1 = 0.0
    time2 = 0.0
    for j = 1:maxJ
        A1 = randn(n, n)
        B1 = randn(n, n)
        A2 = deepcopy(A1)
        B2 = deepcopy(B1)
        ans1 = @timed svd!(A1, B1)
        ans2 = @timed gsvdwrapper(A2, B2)
        time1 = time1 + ans1[2]
        time2 = time2 + ans2[2]
    end
        tjsvd = time1/maxJ
        tjwrapper = time2/maxJ
    return tjsvd, tjwrapper
end
