using Plots
using DelimitedFiles
using .GSVD
using LinearAlgebra

function genetest(m, p, n)
    # y = readdlm("m-p-n=$m-$p-$n.txt")
    maxI = 20
    maxJ = 30
    x = zeros(maxI)
    tsvd = zeros(maxI)
    tgsvd = zeros(maxI)

    for i = 1:maxI
        x[i] = m*i*6
        time1 = 0.0
        time2 = 0.0
        for j = 1:maxJ
            A1 = randn(m*i*6, n*i*6)
            B1 = randn(p*i*6, n*i*6)
            A2 = deepcopy(A1)
            B2 = deepcopy(B1)
            ans1 = @timed svd(A1, B1)
            ans2 = @timed gsvd(A2, B2)
            time1 = time1 + ans1[2]
            time2 = time2 + ans2[2]
        end
            tsvd[i] = time1/maxJ
            tgsvd[i] = time2/maxJ

        println("m: ", x[i])
        println("--------------------------------------------------------")
    end

    y = [tsvd tgsvd];

    display(y);

    x = zeros(maxI,1)
    for i=1:maxI
        x[i] = m*6*i
    end
    gr();
    pl = plot(x, y,
    # xaxis=:log,
    yaxis=:log,
    title = "CPU timing of GSVD, m:p:n = $m:$p:$n",
    label = ["svd(A, B)" "gsvd(A, B)"],
    # yformatter = :scientific,
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
    savefig(pl, "m-p-n=$m-$p-$n.pdf")
    savefig("m-p-n=$m-$p-$n.png")

    open("m-p-n=$m-$p-$n.txt", "w") do io
        writedlm(io, y)
    end

    pl = bar(
    string.(6*m:6*m:6*m*maxI),
    y[:,1]./y[:,2],
    xlabel = "m: number of rows of A",
    xticks = :all,
    ylabel = "factors of speedup",
    title = "t(svd(A, B))/t(gsvd(A, B)), m:p:n=$m:$p:$n",
    fillcolor = :royalblue1,
    fillalpha=0.25,
    legend = nothing)
    savefig(pl, "speedup_m-p-n=$m-$p-$n.pdf")
    savefig("speedup_m-p-n=$m-$p-$n.png")
end
