# ENV["MPLBACKEND"]="qt5agg";
include("../src/gsvd.jl")
# using PyPlot
# pygui(true)
using Plots
using DelimitedFiles

function jmtest()
    maxI = 20
    maxJ = 15
    x = zeros(maxI)
    tjulia = zeros(maxI)

    # pwd() to check current working directory
    tmatlab = readdlm("MATLAB_timing.txt")

    for i = 1:maxI
        x[i] = 100*i;
        println("dimension: ", 100*i);
        t = zeros(maxJ)
        for j = 1:maxJ
            A = randn(100*i, 100*i)
            B = randn(100*i, 100*i)
            ans = @timed gsvd(A, B)
            t[j] = ans[2]
        end
            tjulia[i] = median(t)
    end

    # plot(dims, tjulia, color="red", linewidth=2.0, linestyle="-", marker="o", label="gsvd in Julia (proposed)");
    # plot(dims, tmatlab, color="blue", linewidth=2.0, linestyle="--", marker="*", label="gsvd in MATLAB 2019b");
    # xlabel("# row of A");
    # ylabel("elapsed time in seconds");
    # title("cpu timing of gsvd, m=p=n");
    # legend();
    # show()

    y = [tmatlab tjulia];
    gr();
    pl = plot(x, y,
    title = "CPU timing of GSVD, m=p=n",
    label = ["M-GSVD" "J-GSVD"],
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
    open("J-M-m=p=n.txt", "w") do io
        writedlm(io, y)
    end
    savefig(pl, "J-M-m=p=n.pdf")
    savefig("J-M-m=p=n.png")
end

function barplot()
    y = readdlm("J-M-m=p=n.txt")
    maxI = 20
    gr();
    pl = bar(
    string.(100:100:100*maxI),
    y[:,2]./y[:,1],
    xlabel = "m: number of rows of A",
    xticks = :all,
    xrotation = 60,
    ylabel = "factors of speedup",
    title = "t(J-GSVD)/t(M-GSVD), m=p=n",
    fillcolor = :royalblue1,
    fillalpha=0.25,
    legend = nothing)
    savefig(pl, "speedup_J-M-m=p=n.pdf")
    savefig("speedup_J-M-m=p=n.png")
end
