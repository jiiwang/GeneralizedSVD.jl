include("../src/wrapper.jl")
include("../src/safeDiag.jl")

using DataFrames
# using DelimitedFiles
using CSV

function test(m, l)
    if m > l
        throw(ArgumentError("m <= l is required"))
    else
        F = qr(randn(m+l, l))
        Q = Matrix(F.Q)
        # cannot have @views
        X11 = Q[1:m,:]
        X21 = Q[m+1:end,:]
        A = deepcopy(X11)
        B = deepcopy(X21)

        @time theta, U1, U2, V1T = dorcsd2by1!(X11, X21)
        r = min(m, l)
        if m == l
            C1 = Diagonal(cos.(theta[1:r]))
            S1 = Diagonal(sin.(theta[1:r]))
        else
            C1 = [Diagonal(cos.(theta[1:r])) zeros(m, l-m)]
            S1 = Matrix(1.0I, l, l)
            @views S1[1:r,1:r] = Diagonal(sin.(theta[1:r]))
        end

        # return U1, U2, V1T, A, B, C1, S1

        e = eps(Float64)
        res_x11 = opnorm(U1 * C1 * V1T - A, 1)/(l*opnorm(A, 1)*e)
        res_x21 = opnorm(U2 * S1 * V1T - B, 1)/(l*opnorm(B, 1)*e)
        orthog_cs = opnorm(C1' * C1 + S1' * S1 - I, 1)/(l*e)
        orthog_u1 = opnorm(U1' * U1 - I, 1)/(m*e)
        orthog_u2 = opnorm(U2' * U2 - I, 1)/(l*e)
        orthog_v1 = opnorm(V1T * V1T' - I, 1)/(l*e)
        return res_x11, res_x21, orthog_cs, orthog_u1, orthog_u2, orthog_v1
    end
end

function test()
    df = DataFrame(m = Int64[], l = Int64[], res_x11 = Float64[],
    res_x21 = Float64[], orthog_cs = Float64[], orthog_u1 = Float64[],
    orthog_u2 = Float64[], orthog_v1 = Float64[])

    # case 1: m < l
    println("===============================================================")
    for i = 1:20
        a, b, c, d, e, f = test(40, 60)
        push!(df, [40 60 a b c d e f])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f = test(200, 300)
        push!(df, [200 300 a b c d e f])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f = test(600, 900)
        push!(df, [600 900 a b c d e f])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f = test(1000, 1500)
        push!(df, [1000 1500 a b c d e f])
        println("---------------------------------------------------------------")
    end
    #
    # case 2: m = l
    println("===============================================================")
    for i = 1:20
        a, b, c, d, e, f = test(50, 50)
        push!(df, [50 50 a b c d e f])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f = test(250, 250)
        push!(df, [250 250 a b c d e f])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f = test(750, 750)
        push!(df, [750 750 a b c d e f])
        println("---------------------------------------------------------------")
    end
    for i = 1:20
        a, b, c, d, e, f = test(1250, 1250)
        push!(df, [1250 1250 a b c d e f])
        println("---------------------------------------------------------------")
    end

    CSV.write("rand_matrix_test_lapack.csv", df, append=true)
    return df
end

# Example 2 in Chapter 3 of Jenny's thesis
function test2()
    Q1 = [1/sqrt(7) 0 0 1/sqrt(3);
          1/sqrt(7) -2/sqrt(10) -1/2 -1/(2*sqrt(3));
          1/sqrt(7) -1/sqrt(10) 1/4  3/(4*sqrt(3))]
    Q2 = [1/sqrt(7) 0 3/4 -3/(4*sqrt(3));
          1/sqrt(7) 0 0 0;
          1/sqrt(7) 1/sqrt(10) -1/4 -3/(4*sqrt(3));
          1/sqrt(7) 2/sqrt(10) -1/4 1/(4*sqrt(3))]
    m = 3; l = 4;
    A = deepcopy(Q1); B = deepcopy(Q2);
    theta, U1, U2, V1T = dorcsd2by1!(Q1, Q2)
    r = min(m, l)
    if m == l
        C1 = Diagonal(cos.(theta[1:r]))
        S1 = Diagonal(sin.(theta[1:r]))
    else
        C1 = [Diagonal(cos.(theta[1:r])) zeros(m, l-m)]
        S1 = Matrix(1.0I, l, l)
        @views S1[1:r,1:r] = Diagonal(sin.(theta[1:r]))
    end
    e = eps(Float64)
    res_x11 = opnorm(U1 * C1 * V1T - A, 1)/(l*opnorm(A, 1)*e)
    res_x21 = opnorm(U2 * S1 * V1T - B, 1)/(l*opnorm(B, 1)*e)
    orthog_cs = opnorm(C1' * C1 + S1' * S1 - I, 1)/(l*e)
    orthog_u1 = opnorm(U1' * U1 - I, 1)/(m*e)
    orthog_u2 = opnorm(U2' * U2 - I, 1)/(l*e)
    orthog_v1 = opnorm(V1T * V1T' - I, 1)/(l*e)
    println("res_q1: ", res_x11)
    println("res_q2: ", res_x21)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u1)
    println("orthog_v: ", orthog_u2)
    println("orthog_z: ", orthog_v1)
    return U1, U2, V1T', C1, S1
end
