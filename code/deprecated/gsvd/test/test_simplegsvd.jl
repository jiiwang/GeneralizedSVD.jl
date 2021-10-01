# ENV["MPLBACKEND"]="qt5agg";
# using PyPlot
# using BenchmarkTools
# pygui(true)
include("../src/simple_gsvd.jl")

using Plots
using DelimitedFiles

function test1()
    A = [1. 4 1 0;5 3 1 1;3 0 1 2;7 2 5 1];
    B = [4. 5 1 3;-2 0 1 4;3 2 1 -5;1 1 -6 3];
    n = size(A)[1]
    A_ = copy(A);
    B_ = copy(B);
    U,V,X,C,S = simplegsvd(A,B);
    lambda = diag(C'*C)./diag(S'*S);
    println(sqrt.(lambda))
    X_ = inv(X');
    for i = 1:n
        # println(A_'*A_*X_[:,i] - lambda[i]*B_'*B_*X_[:,i])
        println(norm(A_'*A_*X_[:,i] - lambda[i]*B_'*B_*X_[:,i],2))
    end
end

function test2()
    A = [1. 2 3 0 6; 5 4 2 1 5; 0 3 5 2 2; 2 1 3 3 4; 2 0 5 3 1];
    B = [1. 0 3 -1 2; -2 5 0 1 1; 4 2 -1 2 4;1 3 0 4 2;-3 5 7 2 8];
    n = size(A)[1]
    A_ = copy(A);
    B_ = copy(B);
    U,V,X,C,S = simplegsvd(A,B);
    lambda = diag(C'*C)./diag(S'*S);
    println(sqrt.(lambda))
    X_ = inv(X');
    for i = 1:n
        # println(A_'*A_*X_[:,i] - lambda[i]*B_'*B_*X_[:,i])
        println(norm(A_'*A_*X_[:,i] - lambda[i]*B_'*B_*X_[:,i],2))
    end
end

function test3()
    A = [0.  1  2  2; 1  2  1  2; 1  1  1  2; 1  2  2  3];
    B = [0.  0  1  1; 1  2  2  3; 1  1  1  2; 1  1  1  2];
    n = size(A)[1]
    A_ = copy(A);
    B_ = copy(B);
    U,V,X,C,S = simplegsvd(A,B);
    println(sqrt.(diag(C'*C)./diag(S'*S)));
    e = eps(Float64)
    res_a = opnorm(A_ - U*C*X', 1)/(n*opnorm(A_, 1)*e)
    res_b = opnorm(B_ - V*S*X', 1)/(n*opnorm(B_, 1)*e)
    orthog_cs = opnorm(I - C'*C - S'*S, 1)/(n*e)
    orthog_u = opnorm(I - U'*U, 1)/(n*e)
    orthog_v = opnorm(I - V'*V, 1)/(n*e)
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    return U,V,X,C,S
end

function test4()
    A =  [4.  2  2  0  0; 0  0  2  0  0; 4  4  4  0  2; 2  2  2  0  0; 2  4  4  0  2];
    B =  [4.  4  2  0  2; 2  2  0  0  2; 2  2  0  0  2; 4  2  2  0  0; 2  0  0  0  0];
    n = size(A)[1]
    A_ = copy(A);
    B_ = copy(B);
    U,V,X,C,S = simplegsvd(A,B);
    println(sqrt.(diag(C'*C)./diag(S'*S)));
    e = eps(Float64)
    res_a = opnorm(A_ - U*C*X', 1)/(n*opnorm(A_, 1)*e)
    res_b = opnorm(B_ - V*S*X', 1)/(n*opnorm(B_, 1)*e)
    orthog_cs = opnorm(I - C'*C - S'*S, 1)/(n*e)
    orthog_u = opnorm(I - U'*U, 1)/(n*e)
    orthog_v = opnorm(I - V'*V, 1)/(n*e)
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    return U,V,X,C,S
end

function jmtest()
    maxI = 20
    maxJ = 15
    x = zeros(maxI)

    tmatlab = readdlm("MATLAB_timing.txt")
    tjulia = zeros(maxI)
    # tmatlab = [0.0027, 0.0087, 0.0198, 0.0341, 0.0552,
    #            0.0809, 0.1156, 0.1571, 0.2225, 0.3001,
    #            0.3884, 0.4945, 0.6240, 0.8282, 1.0740,
    #            1.2793, 2.4995, 2.1374, 2.5280, 3.1289]

    # tmatlab = [0.0024, 0.0086, 0.0196, 0.0352, 0.0561,
    #            0.0845, 0.1171, 0.1610, 0.2244, 0.2936,
    #            0.3875, 0.5220, 0.6757, 0.8499, 1.1412,
    #            1.5362, 2.0623, 2.1121, 2.5934, 3.1820]

    for i = 1:maxI
        x[i] = 100*i;
        println("dimension: ", 100*i)
        t = zeros(maxJ)
        for j = 1:maxJ
            A = randn(100*i, 100*i)
            B = randn(100*i, 100*i)
            ans = @timed simplegsvd(A, B)
            t[j] = ans[2]
        end
        tjulia[i] = median(t)
        # A = randn(100*i, 100*i)
        # B = randn(100*i, 100*i)
        # bench = @benchmark simplegsvd($A, $B)
        # tjulia[i] = mean(bench.times)
    end

    # print(tjulia)
    # [0.0038232673999999997, 0.015201680799999999, 0.02848866366666666, 0.0540633, 0.08265005580000001,
    # 0.1193703124, 0.16424617279999998, 0.2354030057333333, 0.3361372554666667, 0.4108958160666667,
    # 0.5092602177333333, 0.6309922593333332, 0.7693808228666666, 0.9428725853333333, 1.3132644359333334,
    # 1.9010355676000001, 2.0724060565999998, 2.540632365533333, 2.9474200295999995, 3.537733264866666]

    # plot(dims, tjulia, color="red", linewidth=2.0, linestyle="-", marker="o", label="simplified-gsvd in Julia");
    # plot(dims, tmatlab, color="blue", linewidth=2.0, linestyle="--", marker="*", label="gsvd in MATLAB 2019b");
    # xlabel("# row of A");
    # ylabel("elapsed time in seconds");
    # title("cpu timing of gsvd, m=p=n");
    # legend();
    # show()
      y = [tjulia tmatlab];
      gr();
      pl = plot(x, y,
      title = "CPU timing of GSVD, m=p=n",
      label = ["SJ-GSVD" "M-GSVD"],
      xlabel = "m: number of rows of A",
      ylabel = "elapsed time in seconds",
      linecolor = [:firebrick2 :steelblue2],
      markershape = [:circle :xcross],
      markercolor = [:firebrick2 :steelblue2],
      markersize = 3.5,
      markeralpha = 0.75,
      markerstrokewidth=0,
      lw = 1.5,
      legend = :topleft)
      # savefig("myplot.png") # Saves the CURRENT_PLOT as a .png
      # savefig(pl, "myplot.pdf") # Saves the plot from p as a .pdf vector graphic
      # return y, pl
      # Check pwd()
      open("SJ-M-m=p=n.txt", "w") do io
          writedlm(io, y)
      end
      savefig(pl, "SJ-M-m=p=n.pdf")
      savefig("SJ-M-m=p=n.png")
end

function profiling()
    for i = 5:5
        tqr = 0
        tq = 0
        tsvd = 0
        tsafediag = 0
        tall = 0
        for j = 1:10
            A = randn(100*i, 100*i)
            B = randn(100*i, 100*i)
            ans = @timed simplegsvd(A, B)
            tqr = tqr + ans[1][6]
            tq = tq + ans[1][7]
            tsvd = tsvd + ans[1][8]
            tsafediag = tsafediag + ans[1][9]
            tall = tall + ans[2]
        end
        tqr = tqr/10
        tq = tq/10
        tsvd = tsvd/10
        tsafediag = tsafediag/10
        tall = tall/10
        println(100*i, " & ", round(tqr, digits=4), " & ", string(round(tqr/tall*100, digits=2),"%"), " & ",
        round(tq, digits=4), " & ", string(round(tq/tall*100, digits=2),"%"), " & ",
        round(tsvd, digits=4), " & ", string(round(tsvd/tall*100, digits=2),"%"), " & ",
        round(tsafediag, digits=4)," & ", string(round(tsafediag/tall*100, digits=2),"%"), " & ",
        string(round(tsvd/tsafediag*100, digits=2),"%"), " & ",
        round(tall, digits=4))
    end
end


# julia> test1()
# [-1.3589129821411916e-13, -8.482103908136196e-14, 5.329070518200751e-15, -7.394085344003543e-14]
# 1.765126125838869e-13
# [8.881784197001252e-16, -2.220446049250313e-15, 1.5543122344752192e-15, -1.7208456881689926e-15]
# 3.331131634682349e-15
# [-6.661338147750939e-16, 3.885780586188048e-16, -7.771561172376096e-16, -1.0547118733938987e-15]
# 1.5202354861220295e-15
# [-4.163336342344337e-16, -6.869504964868156e-16, -1.096345236817342e-15, 3.0531133177191805e-16]
# 1.3929905122500834e-15
#
# julia> test2()
# [1.1546319456101628e-14, 4.1744385725905886e-14, 1.509903313490213e-14, 2.1316282072803006e-14, 7.105427357601002e-14]
# 8.721803545738815e-14
# [-2.220446049250313e-15, 1.3322676295501878e-15, -4.218847493575595e-15, 1.6653345369377348e-15, -4.440892098500626e-15]
# 6.855570991454365e-15
# [-3.1086244689504383e-15, 0.0, 1.7763568394002505e-15, 4.440892098500626e-16, -3.552713678800501e-15]
# 5.063396036227354e-15
# [1.9984014443252818e-15, 9.992007221626409e-16, -2.6645352591003757e-15, -8.465450562766819e-16, -4.440892098500626e-16]
# 3.606329480434721e-15
# [5.828670879282072e-16, -1.1102230246251565e-16, 2.2065682614424986e-15, 1.0269562977782698e-15, 1.0547118733938987e-15]
# 2.718098574309599e-15
