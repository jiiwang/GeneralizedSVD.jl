ENV["MPLBACKEND"]="qt5agg";
using PyPlot;
pygui(true)

function genetest(m, p, n)
    # params:
    # m : row(A)
    # p : row(B)
    # n : col(A) or col(B)

    # maxI: # of tests
    maxI = 15
    # maxJ: # of runs
    maxJ = 10
    X = zeros(maxI)
    tgsvd = zeros(maxI)
    tsvd = zeros(maxI)

    for i = 1:5
        X[i] = 10*m*i*2
        time1 = 0.0
        time2 = 0.0
        for j = 1:maxJ
            A1 = randn(10*m*i*2, 10*n*i*2)
            B1 = randn(10*p*i*2, 10*n*i*2)
            A2 = deepcopy(A1)
            B2 = deepcopy(B1)
            ans1 = @timed gsvd(A1, B1)
            ans2 = @timed svd(A2, B2)
            time1 = time1 + ans1[2]
            time2 = time2 + ans2[2]
        end
            tgsvd[i] = time1/maxJ
            tsvd[i] = time2/maxJ

        # e = eps(Float64)
        # res_a1 = opnorm(ans1[1][1]'*A_*ans1[1][3] - ans1[1][4]*ans1[1][6], 1)/(max(10*m*i,10*n*i)*opnorm(A_, 1)*e)
        # res_a2 = opnorm(ans2[1].U'*A_*ans2[1].Q - ans2[1].D1*ans2[1].R0, 1)/(max(10*m*i,10*n*i)*opnorm(A_, 1)*e)
        println("m: ", X[i])
        # println("residual err of A by new gsvd: ", res_a1)
        # println("residual err of A by current gsvd in Julia: ", res_a2)
        println("--------------------------------------------------------")
    end

    for i = 6:maxI
        X[i] = 10*m*i*2
        time1 = 0.0
        time2 = 0.0
        for j = 1:2
            A1 = randn(10*m*i*2, 10*n*i*2)
            B1 = randn(10*p*i*2, 10*n*i*2)
            A2 = deepcopy(A1)
            B2 = deepcopy(B1)
            ans1 = @timed gsvd(A1, B1)
            ans2 = @timed svd(A2, B2)
            time1 = time1 + ans1[2]
            time2 = time2 + ans2[2]
        end
            tgsvd[i] = time1/maxJ
            tsvd[i] = time2/maxJ

        # e = eps(Float64)
        # res_a1 = opnorm(ans1[1][1]'*A_*ans1[1][3] - ans1[1][4]*ans1[1][6], 1)/(max(10*m*i,10*n*i)*opnorm(A_, 1)*e)
        # res_a2 = opnorm(ans2[1].U'*A_*ans2[1].Q - ans2[1].D1*ans2[1].R0, 1)/(max(10*m*i,10*n*i)*opnorm(A_, 1)*e)
        println("m: ", X[i])
        # println("residual err of A by new gsvd: ", res_a1)
        # println("residual err of A by current gsvd in Julia: ", res_a2)
        println("--------------------------------------------------------")
    end

    println("native gsvd(proposed)")
    println(tgsvd)
    println("gsvd in Julia 1.3")
    println("--------------------------------------------------------")
    println(tsvd)

    plot(X, tgsvd, color="red", linewidth=2.0, linestyle="-", marker="o",label="native gsvd (proposed)");
    plot(X, tsvd, color="blue", linewidth=2.0, linestyle="--", marker="*", label="gsvd in Julia 1.3");
    xlabel("# row of A");
    ylabel("elapsed time in seconds");
    title("cpu timing of gsvd, m:p:n = $m:$p:$n");
    legend();
    show()
end
