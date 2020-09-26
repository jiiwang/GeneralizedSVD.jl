using DelimitedFiles

function metrictest()
    cd("/Users/hytonwons/Study/Julia-Workspace/GSVD_julia")
    return readfile()
    # generator()
    # writefile()
end

function readfile()
    F = readdlm("gsvd/metric_dimension.csv", ',', Int, '\n', skipstart=1)
end
