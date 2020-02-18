ENV["MPLBACKEND"]="qt5agg";
using PyPlot;
pygui(true)

function genetest(m, n, p)
    # params: row(A): row(B): col ratio

    maxI = 20
    X = zeros(maxI)
    tgsvd = zeros(maxI)
    tsvd = zeros(maxI)

    for i = 1:maxI
        X[i] = 10*m*i
        A1 = randn(10*m*i, 10*n*i)
        B1 = randn(10*p*i, 10*n*i)
        A2 = copy(A1)
        B2 = copy(B1)
        A_ = copy(A1)
        B_ = copy(B1)
        ans1 = @timed gsvd(A1, B1, 1)
        ans2 = @timed svd(A2, B2)
        tgsvd[i] = ans1[2]
        tsvd[i] = ans2[2]

        e = eps(Float64)
        res_a1 = norm(ans1[1][1]'*A_*ans1[1][3] - ans1[1][4]*ans1[1][6], 1)/(max(10*m*i,10*n*i)*norm(A_, 1)*e)
        res_a2 = norm(ans2[1].U'*A_*ans2[1].Q - ans2[1].D1*ans2[1].R0, 1)/(max(10*m*i,10*n*i)*norm(A_, 1)*e)
        println("m: ", X[i])
        println("residual err of A by new gsvd: ", res_a1)
        println("residual err of A by current gsvd in Julia: ", res_a2)
        println("--------------------------------------------------------")
    end

    plot(X, tgsvd, color="red", linewidth=2.0, linestyle="-", label="new gsvd");
    plot(X, tsvd, color="blue", linewidth=2.0, linestyle="--", label="current gsvd");
    xlabel("# row of A");
    ylabel("elapsed time in seconds");
    title("cpu timing of gsvd, m:n:p = $m:$n:$p");
    legend();
    show()
end
