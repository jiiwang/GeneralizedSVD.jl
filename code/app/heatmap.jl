using Plots
include("../gsvd/src/gsvd.jl")

Y = readdlm("Y.txt");
H = readdlm("H.txt");
# F = gsvd(Y, H);

colGRAD = cgrad([colorant"green",colorant"black",colorant"red"]);

function heatmapYeast()
    # xs = [string(i) for i = 1:18];
    # contour(Y, fill=true, lw=0, c=colGRAD, clim=(0,3))
    # heatmap(xs, Y, c=colGRAD, clim=(0,3))
    heatmap([string(i) for i = 1:18], [string(i) for i = 1:20], randn(20, 18), yflip=true, c=colGRAD, clim=(0, 1))
end

function heatmapHuman()
end

function heatmapU()
end

function heatmapV()
end
