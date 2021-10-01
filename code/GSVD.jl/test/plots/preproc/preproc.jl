include("../../wrapper/preprocWrapper.jl")
include("../../../src/preproc.jl")

# using BenchmarkTools
using Plots
using DelimitedFiles

function preprocplot(m, p, n)
    # params: row(A): row(B): col ratio

    maxI = 10
    maxJ = 100
    x = zeros(maxI)
    twrapper = zeros(maxI)
    tnative = zeros(maxI)

    for i = 1:maxI
        x[i] = m*i*2
        time1 = 0.0
        time2 = 0.0
        for j = 1:maxJ
            A1 = randn(m*i*2, n*i*2)
            B1 = randn(p*i*2, n*i*2)
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
    yformatter = :scientific,
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

function preprocbarplot(m, p, n)
    # y = randn(20);
    y = readdlm("preproc-m-p-n=$m-$p-$n.txt")
    maxI = 10
    gr();
    pl = bar(
    string.(m*2:m*2:m*maxI*2),
    y[:,1]./y[:,2],
    xlabel = "m: number of rows of A",
    xticks = :all,
    ylabel = "factors of speedup",
    title = "preproc: t(LAPACK-wrapper)/t(Native-Julia), m:p:n=$m:$n:$p",
    fillcolor = :royalblue1,
    fillalpha=0.25,
    legend = nothing)
    savefig(pl, "speedup_preproc-m-p-n=$m-$p-$n.pdf")
    savefig("speedup_preproc-m-p-n=$m-$p-$n.png")
end
