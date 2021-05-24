ENV["MPLBACKEND"]="qt5agg";
include("../gsvd/src/gsvd.jl")
include("../gsvd/test/test.jl")
using PyPlot
pygui(true)
# using Plots
using DataFrames
using DelimitedFiles
using LinearAlgebra:svd
using Statistics:mean

function comparativeAnalyser()
    yeast, human = checkNull()
    A = SVDimputation(yeast)
    B = SVDimputation(human)

    A_ = deepcopy(A)
    B_ = deepcopy(B)
    F = gsvd(A, B)

    res_a, res_b, orthog_cs, orthog_u, orthog_v, orthog_q = computemetric(F, A_, B_)
    println("res_a: ", res_a)
    println("res_b: ", res_b)
    println("orthog_cs: ", orthog_cs)
    println("orthog_u: ", orthog_u)
    println("orthog_v: ", orthog_v)
    println("orthog_q: ", orthog_q)

    dist = atan.(sqrt.(diag(F.D1'*F.D1)./diag(F.D2'*F.D2))) .- pi/4
    println(dist);
    y = ["18","17","16","15","14","13","12","11","10","9","8","7","6","5","4","3","2","1"];
    ax = gca();
    PyPlot.xlim([-pi/4,pi/4]);
    ax.set_xticks(range(-pi/4, pi/4, step=pi/8));
    ax.set_xticklabels([L"-$\frac{\pi}{4}$", L"-$\frac{\pi}{8}$", L"$0$", L"$\frac{\pi}{8}$", L"$\frac{\pi}{4}$"]);
    ax.xaxis.grid(true);
    b = barh(y,dist[end:-1:1],color="#0f87bf",align="center",alpha=0.4);PyPlot.title("Angular Distance")
end

# function matrixRetrieval()
#     yeast, human = checkNull()
#     A = SVDimputation(yeast)
#     B = SVDimputation(human)
#     return A, B
# end


function slicingData()
    human = readdlm("input_data/Human.txt")[3:end,6:end]
    yeast = readdlm("input_data/Yeast.txt")[3:end,7:end]

    open("input_data/human_gene_exp.txt", "w") do io
        writedlm(io, human)
    end

    open("input_data/yeast_gene_exp.txt", "w") do io
        writedlm(io, yeast)
    end
end

# function estimate()
#     H = readdlm("input_data/human_nomissing.txt")
#     F = svd(H, full = true)
#     n = size(H)[2]
#
#     fractions = (F.S.^2) ./ sum(F.S.^2)
#     println(fractions[2])
#     entropy = -(1/log(n)) * sum(fractions .* log.(fractions))
#     println(entropy)
# end

function checkNull()
    yeast = readdlm("input_data/yeast_gene_exp.txt")
    for i=1:size(yeast)[1]
        for j=1:size(yeast)[2]
            if yeast[i,j] == "Null"
                yeast[i,j] = missing
            end
        end
    end

    human = readdlm("input_data/human_gene_exp.txt")
    # count = 0
    for i=1:size(human)[1]
        # if "Null" in human[i,:]
        #     count = count+1
        for j=1:size(human)[2]
            if human[i,j] == "Null"
                human[i,j] = missing
            end
        end
        # end
    end
    return yeast, human
end

function SVDimputation(input)
    # yeast = readdlm("input_data/yeast_readable.txt")

    inputraw = deepcopy(input)

    output = zeros(size(input)[1], size(input)[2])

    # calculate column mean and assign it to the missing entries
    for j=1:size(input)[2]
        m = mean(skipmissing(input[:,j]))
        for i=1:size(input)[1]
            if typeof(input[i,j]) == Missing
                output[i,j] = m
            else
                output[i,j] = convert(Float64, input[i,j])
            end
        end
    end

    maxItr = 10
    tol = 1.0e-10

    i = 1
    sig = 5
    err = 1.0
    while i <= maxItr && err > tol
        F = svd(output, full = true)
        approx = F.U[:,1:sig]*Diagonal(F.S[1:sig])*F.Vt[1:sig,:]
        outputnew = deepcopy(output)
        for j=1:size(inputraw)[2]
            for i=1:size(inputraw)[1]
                if typeof(inputraw[i,j]) == Missing
                    outputnew[i,j] = approx[i,j]
                end
            end
        end
        err = opnorm(output - outputnew, 1)/opnorm(output, 1)
        output = deepcopy(outputnew)
        i = i+1
    end
    println("# iterations: ", i)
    return output
end
