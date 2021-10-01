include("generatebyrank.jl")
include("../src/preproc.jl")
include("../src/wrapper.jl")

using Plots
using DelimitedFiles

function test(m, p, n)
    A = randn(m, n)
    B = randn(p, n)
    A1 = deepcopy(A)
    B1 = deepcopy(B)
end

function test1()
    A1 = [1.0 2 1 0;2 3 1 1;3 4 1 2;4 5 1 3;5 6 1 4]
    B1 = [6.0 7 1 5;7 1 -6 13;-4 8 9 -2]

    A2 = [1.0 2 1 0;2 3 1 1;3 4 1 2;4 5 1 3;5 6 1 4]
    B2 = [6.0 7 1 5;7 1 -6 13;-4 8 9 -2]

    A = [1.0 2 1 0;2 3 1 1;3 4 1 2;4 5 1 3;5 6 1 4]
    B = [6.0 7 1 5;7 1 -6 13;-4 8 9 -2]
    m, n = size(A1)
    p = size(B1)[1]

    U1, V1, Q1, k1, l1, A1, B1 = @time dggsvp3!(A1, B1)
    println("k: ", k1)
    println("l: ", l1)
    U2, V2, Q2, k2, l2, A2, B2 = @time preproc(A2, B2)
    # A2, U, Q = preproc(A2, B2)

    e = eps(Float64)
    println("Stability test on dggsvp3")
    println("res_a = ", opnorm(U1*A1*Q1' - A, 1)/(max(m,n)*opnorm(A, 1)*e))
    println("res_b = ", opnorm(V1*B1*Q1' - B, 1)/(max(m,n)*opnorm(B, 1)*e))
    println("orthog_u = ", opnorm(U1'*U1-I, 1)/(m*e))
    println("orthog_v = ", opnorm(V1'*V1-I, 1)/(p*e))
    println("orthog_q = ", opnorm(Q1'*Q1-I, 1)/(n*e))
    println("--------------------------------------")
    println("Stability test on preproc")
    println("res_a = ", opnorm(U2*A2*Q2' - A, 1)/(max(m,n)*opnorm(A, 1)*e))
    println("res_b = ", opnorm(V2*B2*Q2' - B, 1)/(max(m,n)*opnorm(B, 1)*e))
    println("orthog_u = ", opnorm(U2'*U2-I, 1)/(m*e))
    println("orthog_v = ", opnorm(V2'*V2-I, 1)/(p*e))
    println("orthog_q = ", opnorm(Q2'*Q2-I, 1)/(n*e))

end

function test2()
    A = [1.0 2 1 0;2 3 1 1;3 4 1 2;4 5 1 3;5 6 1 4]
    B = [6.0 7 1 5;7 1 -6 13;-4 8 9 -2]
    return preproc(A,B)
end

function test3(m, p, n)
    A = randn(m, n)
    B = randn(p, n)
    A1 = deepcopy(A)
    B1 = deepcopy(B)
    A2 = deepcopy(A)
    B2 = deepcopy(B)
    U1, V1, Q1, k1, l1, A1, B1 = @time dggsvp3!(A1, B1)
    println("k: ", k1)
    println("l: ", l1)
    U2, V2, Q2, k2, l2, A2, B2 = @time preproc(A2, B2)
    # A2, U, Q = preproc(A2, B2)

    e = eps(Float64)
    println("Stability test on dggsvp3")
    println("res_a = ", opnorm(U1*A1*Q1' - A, 1)/(max(m,n)*opnorm(A, 1)*e))
    println("res_b = ", opnorm(V1*B1*Q1' - B, 1)/(max(m,n)*opnorm(B, 1)*e))
    println("orthog_u = ", opnorm(U1'*U1-I, 1)/(m*e))
    println("orthog_v = ", opnorm(V1'*V1-I, 1)/(p*e))
    println("orthog_q = ", opnorm(Q1'*Q1-I, 1)/(n*e))
    println("--------------------------------------")
    println("Stability test on preproc")
    println("res_a = ", opnorm(U2*A2*Q2' - A, 1)/(max(m,n)*opnorm(A, 1)*e))
    println("res_b = ", opnorm(V2*B2*Q2' - B, 1)/(max(m,n)*opnorm(B, 1)*e))
    println("orthog_u = ", opnorm(U2'*U2-I, 1)/(m*e))
    println("orthog_v = ", opnorm(V2'*V2-I, 1)/(p*e))
    println("orthog_q = ", opnorm(Q2'*Q2-I, 1)/(n*e))

    return k2, l2
end


function preprocplot(m, p, n)
    # params: row(A): row(B): col ratio

    maxI = 15
    maxJ = 20
    x = zeros(maxI)
    twrapper = zeros(maxI)
    tnative = zeros(maxI)

    for i = 1:maxI
        x[i] = 10*m*i*2
        time1 = 0.0
        time2 = 0.0
        for j = 1:maxJ
            A1 = randn(10*m*i*2, 10*n*i*2)
            B1 = randn(10*p*i*2, 10*n*i*2)
            A2 = deepcopy(A1)
            B2 = deepcopy(B1)
            ans1 = @timed dggsvp3!(A1, B1)
            ans2 = @timed preproc!(A2, B2)
            time1 = time1 + ans1[2]
            time2 = time2 + ans2[2]
        end
            twrapper[i] = time1/maxJ
            tnative[i] = time2/maxJ

        # e = eps(Float64)
        # res_a1 = opnorm(ans1[1][1]'*A_*ans1[1][3] - ans1[1][4]*ans1[1][6], 1)/(max(10*m*i,10*n*i)*opnorm(A_, 1)*e)
        # res_a2 = opnorm(ans2[1].U'*A_*ans2[1].Q - ans2[1].D1*ans2[1].R0, 1)/(max(10*m*i,10*n*i)*opnorm(A_, 1)*e)
        println("m: ", x[i])
        # println("residual err of A by new gsvd: ", res_a1)
        # println("residual err of A by current gsvd in Julia: ", res_a2)
        println("--------------------------------------------------------")
    end

    y = [twrapper tnative];
    gr();
    pl = plot(x, y,
    title = "CPU timing of preproc, m:p:n = $m:$p:$n",
    label = ["LAPACK-wrapper" "Native-Julia"],
    xlabel = "m: number of rows of A",
    ylabel = "elapsed time in seconds",
    markershape = [:circle :xcross],
    markersize = 3.5,
    markeralpha = 0.75,
    markerstrokewidth=0,
    lw = 1.5,
    legend = :topleft)
    # savefig("myplot.png") # Saves the CURRENT_PLOT as a .png
    # savefig(pl, "myplot.pdf") # Saves the plot from p as a .pdf vector graphic
    # return y, pl
    open("preproc-m-p-n=$m-$p-$n.txt", "w") do io
        writedlm(io, y)
    end
    savefig(pl, "preproc-m-p-n=$m-$p-$n.pdf")
    savefig("preproc-m-p-n=$m-$p-$n.png")
end

function barplot(m, n, p)
    # y = randn(20);
    y = readdlm("preproc-m-p-n=$m-$p-$n.txt")
    maxI = 15
    gr();
    pl = bar(
    string.(10*m*2:10*m*2:10*m*maxI*2),
    y[:,1]./y[:,2],
    xlabel = "m: number of rows of A",
    xticks = :all,
    ylabel = "factors of speedup",
    title = "preproc: t(Julia 1.3)/t(LAPACK-Julia X.X), m:p:n=$m:$n:$p",
    fillcolor = :royalblue1,
    fillalpha=0.25,
    legend = nothing)
    savefig(pl, "speedup_preproc-m-p-n=$m-$p-$n.pdf")
    savefig("speedup_preproc-m-p-n=$m-$p-$n.png")
end


function demoplot()
    x = 1:10; y = rand(10, 2); # 2 columns means two lines
    # plotly()
    gr();
    plot(x, y,
    title = "CPU timing of preproc, m:p:n",
    label = ["LAPACK-wrapper" "Native-Julia"],
    xlabel = "m: # row of A",
    ylabel = "elapsed time in seconds",
    linecolor = [:firebrick2 :steelblue2],
    markershape = [:circle :xcross],
    markercolor = [:firebrick2 :steelblue2],
    markersize = 3.5,
    markeralpha = 0.75,
    markerstrokewidth=0,
    lw = 1.5,
    legend = :topleft)
end

# test timing of preproc
# m = 300
# n = 300
# p = 200
# k = 80
# l = 200
function test_profiling()
    A11 = generator(300, 100, 80);
    A12 = randn(300, 200);
    A = [A11 A12];

    println("rank of A11: ",rank(A11));

    F = svd(A11);


    tola = max(300,300)*opnorm(A, 1)*eps(Float64);
    tola1 = max(300,300)*eps(opnorm(A, 1));
    println(tola);
    println(tola1);

    println(diag(A)[1:100])

    k = 0;
    for i = 1:min(300, 100)
        if abs(A[i,i]) > tola1
            k += 1
        end
    end
    println("k: ", k)

    k1 = 0;
    for i = 1:size(F.S)
        if S[i] > tola
            k1 += 1
        end
    end

    println("k1: ", k+1)

    B = generator(200, 300, 200);
    println("l: ", rank(B));
    # @time preproc(A,B);
end

# Example 3 in the special case note
function test4()
    A = [0.  1  2  2; 1  2  1  2; 1  1  1  2; 1  2  2  3];
    B = [0.  0  1  1; 1  2  2  3; 1  1  1  2; 1  1  1  2];
    n = size(A)[1]
    U, V, Q, k, l, A1, B1 = preproc(A, B)
    println("Stability test on preproc")
    e = eps(Float64)
    println("res_a = ", opnorm(U*A1*Q' - A, 1)/(n*opnorm(A, 1)*e))
    println("res_b = ", opnorm(V*B1*Q' - B, 1)/(n*opnorm(B, 1)*e))
    println("orthog_u = ", opnorm(U'*U-I, 1)/(n*e))
    println("orthog_v = ", opnorm(V'*V-I, 1)/(n*e))
    println("orthog_q = ", opnorm(Q'*Q-I, 1)/(n*e))
    return U, V, Q, k, l, A1, B1
end

# Example 4 in the special case note
function test5()
    A =  [4.  2  2  0  0; 0  0  2  0  0; 4  4  4  0  2; 2  2  2  0  0; 2  4  4  0  2];
    B =  [4.  4  2  0  2; 2  2  0  0  2; 2  2  0  0  2; 4  2  2  0  0; 2  0  0  0  0];
    n = size(A)[1]
    U, V, Q, k, l, A1, B1 = preproc(A, B)
    println("Stability test on preproc")
    e = eps(Float64)
    println("res_a = ", opnorm(U*A1*Q' - A, 1)/(n*opnorm(A, 1)*e))
    println("res_b = ", opnorm(V*B1*Q' - B, 1)/(n*opnorm(B, 1)*e))
    println("orthog_u = ", opnorm(U'*U-I, 1)/(n*e))
    println("orthog_v = ", opnorm(V'*V-I, 1)/(n*e))
    println("orthog_q = ", opnorm(Q'*Q-I, 1)/(n*e))
    return U, V, Q, k, l, A1, B1
end
