using LinearAlgebra

# include("../../safediag/src/safeDiag.jl")

# A, B have the same dimensions
function simplegsvd(A, B)
    n = size(A)[1]
    F = qr!([A;B])
    # ans1 = @timed qr!([A;B])
    # F = ans1[1]
    Q = Matrix(F.Q)

    # ans2 = @timed Matrix(F.Q)
    # Q = ans2[1]
    Q1 = @views Q[1:n,:]
    Q2 = @views Q[n+1:end,:]

    # return simplecsd(Q1, Q2)
    U,V,Z,C,S = simplecsd(Q1, Q2)
    # ans3 = @timed simplecsd(Q1, Q2)

    return U, V, (F.R)'*Z, C, S
    # return ans3[1][1], ans3[1][2], (F.R)'*ans3[1][3], ans3[1][4], ans3[1][5],
    # ans1[2], ans2[2], ans3[1][6], ans3[2]
end

function simplecsd(Q1, Q2)
    n = size(Q1)[1]
    F1 = svd!(Q1)
    # ans1 = @timed svd!(Q1)
    # F = ans1[1]
    t = 0
    for i=1:n
        if F1.S[i] ≈ 1.0
            t = t+1
        end
    end
    W = Q2*F1.V

    sine = zeros(n)
    @views sine[t+1:end] = sqrt.((1 .- (F1.S[t+1:end]).^2))
    V = zeros(n, n)
    @views V[:,t+1:n] = W[:,t+1:n] * Diagonal(1 ./ sine[t+1:n])
    if t >= 1
        # @views orth = nullspace(V[:,t+1:n]')
        # for i=1:t
        #     @views V[:,i] = normalize(orth[:,i])
        # end
        @views F2 = qr(V[:,t+1:n])
        @views V[:,1:t] = F2.Q[:,n-t+1:n]
    end

    # while t >= 1
    #     println(t)
    #     @views orth = nullspace(V[:,t+1:n]')
    #     @views V[:,t] = normalize(orth[:,1])
    #     t = t-1
    # end

    return F1.U, V, F1.V, Diagonal(F1.S), Diagonal(sine)
    # , ans1[2]
end
