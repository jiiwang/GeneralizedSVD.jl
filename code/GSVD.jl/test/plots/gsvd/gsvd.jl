using Plots
using DelimitedFiles

function genetest(m, p, n)
    y = readdlm("m-p-n=$m-$p-$n.txt")
    maxI = 15
    x = zeros(maxI,1)
    for i=1:maxI
        x[i] = 10*m*2*i
    end
    gr();
    pl = plot(x, y,
    # xaxis=:log, 
    yaxis=:log,
    title = "CPU timing of L-GSVD, m:p:n = $m:$p:$n",
    label = ["svd(A, B)" "gsvd(A, B)"],
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

    pl = bar(
    string.(10*m*2:10*m*2:10*m*maxI*2),
    y[:,1]./y[:,2],
    xlabel = "m: number of rows of A",
    xticks = :all,
    ylabel = "factors of speedup",
    title = "t(svd(A, B))/t(gsvd(A, B)), m:p:n=$m:$n:$p",
    fillcolor = :royalblue1,
    fillalpha=0.25,
    legend = nothing)
    savefig(pl, "speedup_m-p-n=$m-$p-$n.pdf")
    savefig("speedup_m-p-n=$m-$p-$n.png")
end
