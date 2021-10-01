# ENV["MPLBACKEND"]="qt5agg";
include("../src/gsvd.jl")
# using PyPlot;
# pygui(true)
using Plots
using DelimitedFiles
# using LaTeXStrings

function genetest(m, p, n)
    # params:
    # m : row(A)
    # p : row(B)
    # n : col(A) or col(B)

    # maxI: # of tests
    maxI = 15
    # maxJ: # of runs
    maxJ = 15
    x = zeros(maxI)
    tgsvd = zeros(maxI)
    tsvd = zeros(maxI)

    for i = 1:maxI
        x[i] = 10*m*i*2
        println("m: ", x[i])
        time1 = 0.0
        time2 = 0.0
        for j = 1:maxJ
            println("# run: ", j)
            A = randn(10*m*i*2, 10*n*i*2)
            B = randn(10*p*i*2, 10*n*i*2)
            ans2 = @timed svd(A, B)
            ans1 = @timed gsvd(A, B)
            time1 = time1 + ans1[2]
            time2 = time2 + ans2[2]
        end
        tgsvd[i] = time1/maxJ
        tsvd[i] = time2/maxJ

        # e = eps(Float64)
        # res_a1 = opnorm(ans1[1][1]'*A_*ans1[1][3] - ans1[1][4]*ans1[1][6], 1)/(max(10*m*i,10*n*i)*opnorm(A_, 1)*e)
        # res_a2 = opnorm(ans2[1].U'*A_*ans2[1].Q - ans2[1].D1*ans2[1].R0, 1)/(max(10*m*i,10*n*i)*opnorm(A_, 1)*e)

        # println("residual err of A by new gsvd: ", res_a1)
        # println("residual err of A by current gsvd in Julia: ", res_a2)
        println("--------------------------------------------------------")
    end

    # println("native gsvd(proposed)")
    # println(tgsvd)
    # println("--------------------------------------------------------")
    # println("gsvd in Julia 1.3")
    # println(tsvd)
    #
    # plot(X, tgsvd, color="red", linewidth=2.0, linestyle="-", marker="o",label="native gsvd (proposed)");
    # plot(X, tsvd, color="blue", linewidth=2.0, linestyle="--", marker="*", label="gsvd in Julia 1.3");
    # xlabel("# row of A");
    # ylabel("elapsed time in seconds");
    # title("cpu timing of gsvd, m:p:n = $m:$p:$n");
    # legend();
    # show()
    y = [tsvd tgsvd];
    # y = readdlm("m-p-n=$m-$p-$n.txt")
    # x = 10*m*2:10*m*2:10*m*2*maxI
    gr();
    pl = plot(x, y,
    title = "CPU timing of GSVD, m:p:n = $m:$p:$n",
    label = ["L-GSVD" "J-GSVD"],
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
    open("m-p-n=$m-$p-$n.txt", "w") do io
        writedlm(io, y)
    end
    savefig(pl, "m-p-n=$m-$p-$n.pdf")
    savefig("m-p-n=$m-$p-$n.png")
end

function barplot(m, n, p)
    # y = randn(20);
    # y = readdlm("m-p-n=$m-$p-$n.txt")
    # x = zeros(20,1)
    # for i=1:20
    #     x[i] = 10*m*2*i
    # end
    maxI = 15
    gr();
    pl = bar(
    string.(10*m*2:10*m*2:10*m*maxI*2),
    y[:,1]./y[:,2],
    xlabel = "m: number of rows of A",
    xticks = :all,
    ylabel = "factors of speedup",
    title = "t(L-GSVD)/t(J-GSVD), m:p:n=$m:$n:$p",
    fillcolor = :royalblue1,
    fillalpha=0.25,
    legend = nothing)
    savefig(pl, "speedup_m-p-n=$m-$p-$n.pdf")
    savefig("speedup_m-p-n=$m-$p-$n.png")
end
