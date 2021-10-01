include("../src/wrapper.jl")
include("../src/safeDiag.jl")
using Plots
using DelimitedFiles

function genetest(m, l)
    if m > l
        throw(ArgumentError("m <= l is required"))
    else
        # maxI: # of tests
        maxI = 15
        # maxJ: # of runs
        maxJ = 15
        x = zeros(maxI)
        tjulia = zeros(maxI)
        tlapack = zeros(maxI)

        for i = 1:maxI
            x[i] = 10*m*i*2
            println("m: ", x[i])
            time1 = 0.0
            time2 = 0.0
            for j = 1:maxJ
                println("# run: ", j)
                F = qr(randn(10*m*i*2+10*l*i*2, 10*l*i*2))
                Q = Matrix(F.Q)
                # cannot have @views
                Q1 = Q[1:10*m*i*2,:]
                Q2 = Q[10*m*i*2+1:end,:]
                Q11 = deepcopy(Q1)
                Q12 = deepcopy(Q2)
                ans2 = @timed safeDiag(Q1, Q2)
                ans1 = @timed dorcsd2by1!(Q11, Q12)
                time1 = time1 + ans1[2]
                time2 = time2 + ans2[2]
            end
            tlapack[i] = time1/maxJ
            tjulia[i] = time2/maxJ

            println("--------------------------------------------------------")
        end

        y = [tjulia tlapack];
        # y = readdlm("m-p-n=$m-$p-$n.txt")
        # x = 10*m*2:10*m*2:10*m*2*maxI
        gr();
        pl = plot(x, y,
        title = "CPU timing of CSD, m:l = $m:$l",
        label = ["J-CSD" "L-CSD"],
        xlabel = "m: number of rows of Q1",
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
        open("csd-m-l=$m-$l.txt", "w") do io
            writedlm(io, y)
        end
        savefig(pl, "csd-m-l=$m-$l.pdf")
        savefig("csd-m-l=$m-$l.png")
    end
end

function barplot(m, l)
    y = readdlm("csd-m-l=$m-$l.txt")
    maxI = 15
    gr();
    pl = bar(
    string.(10*m*2:10*m*2:10*m*2*maxI),
    y[:,2]./y[:,1],
    xlabel = "m: number of rows of A",
    xticks = :all,
    xrotation = 60,
    ylabel = "factors of speedup",
    title = "t(J-CSD)/t(L-CSD), m:l=$m:$l",
    fillcolor = :royalblue1,
    fillalpha=0.25,
    legend = nothing)
    savefig(pl, "speedup_csd-m-l=$m-$l.pdf")
    savefig("speedup_csd-m-l=$m-$l.png")
end
