ENV["MPLBACKEND"]="qt5agg";
using PyPlot;
pygui(true)

function jmtest()
    maxI = 20
    maxJ = 10
    dims = zeros(maxI)
    tjulia = zeros(maxI)
    tmatlab = zeros(maxI)
    timeimport = readdlm("MATLAB_timing.csv", ',', Float64)
    for i=1:maxI
        tmatlab[i] = timeimport[i,1]
    end

    for i = 1:maxI
        dims[i] = 100*i;
        println("dimension: ", 100*i)
        for j = 1:maxJ
            A1 = randn(100*i, 100*i)
            B1 = randn(100*i, 100*i)
            ans = @timed gsvd(A1, B1, 1)

            tjulia[i] = tjulia[i] + ans[2]
        end
            tjulia[i] = tjulia[i]/maxJ
    end

    plot(dims, tjulia, color="red", linewidth=2.0, linestyle="-", label="gsvd in Julia");
    plot(dims, tmatlab, color="blue", linewidth=2.0, linestyle="--", label="gsvd in MATLAB");
    xlabel("# row of A");
    ylabel("elapsed time in seconds");
    title("cpu timing of gsvd, m=p=n");
    legend();
    show()
end
