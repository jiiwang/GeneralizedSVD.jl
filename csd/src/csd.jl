# using Pkg
# Pkg.add("LinearAlgebra")
using LinearAlgebra

# This function computes the cosine-sine decomposition (CSD) of an
# (m+p) by l orthogonal matrix Q = (Q1' Q2')', such that:
#
# Q1 = U * C * Z',  Q2 = V * S * Z'
#
# where U, V and Z are orthogonal matrices of size m-by-m, p-by-p,
# and l-by-l, respectively.
# C and S are  m-by-l and p-by-l "diagonal" matrices and
# C' * C + S' * S = I.
#
# C and S have the following structures(both alpha and beta are of length l):
#
# (1) If m >= l and p >= l:
#         C =                           l
#                   l ( diag(alpha(1), ... , alpha(l)) )
#                 m-l (      0                         )
#
#         S =                           l
#                   l ( diag(beta(1), ... , beta(l))   )
#                 p-l (     0                          )
#
# (2) If m >= l and p < l:
#         C =         l-p                        p
#                l-p ( I                         0              )
#                  p ( 0     diag(alpha(l-p+1), ... , alpha(l)) )
#                m-l ( 0                         0              )
#
#         S =         l-p                        p
#                  p ( 0       diag(beta(l-p+1), ... , beta(l)) )
#
# (3) If m <= l and p >= l,
#         C =                            m                 l-m
#                    m ( diag(alpha(1), ... , alpha(m))     0 )
#
#         S =                            m                 l-m
#                   m ( diag(beta(1), ... , beta(m)))       0 )
#                 l-m (                 0                   I )
#                 p-l (                 0                   0 )
#
# (4) If m <= l and p < l,
#          C =            l-p                       m+p-l               l-m
#                   l-p ( I                          0                   0 )
#                 m+p-l ( 0      diag(alpha(l-p+1), ... , alpha(m) )     0 )
#
#          S =            l-p                       m+p-l               l-m
#                 m+p-l ( 0      diag(beta(l-p+1), ... , beta(m) )      0 )
#                   l-p ( 0                         0                   I )
#
# argument
# Q1: m by l matrix
# Q2: p by l matrix
# option: 0 or 1
#
# returns U, V, Z, alpha, beta if option = 0
# U, V, Z, C, S if option = 1

function csd(Q1, Q2, option)
# function csd(Q1::Array{Float64,2}, Q2::Array{Float64,2})
    m = size(Q1)[1] # num of rows in Q1
    p = size(Q2)[1] # num of rows in Q2
    l = size(Q1)[2] # num of cols in Q1/Q2

    # println("# of rows in Q1: ", m)
    # println("# of rows in Q2: ", p)
    # println("# of cols in Q1/Q2: ", l)

    q1 = min(m, l)
    q2 = min(p, l)

    # case 1: m<=p
    if m <= p
        # compute the full SVD of Q2
        D1 = svd(Q2, full = true)

        # swap V's col
        D1.U[:,1:q2] = D1.U[:, q2:-1:1]
        # for k = 1:size(D1.U,1), i = 1:round(Int, q2/2)
        #     D1.U[k,i], D1.U[k,q2+1-i] = D1.U[k,q2+1-i], D1.U[k,i]
        # end
        # swap Z's col or Zt's row
        D1.V[:,1:l] = D1.V[:,l:-1:1]
        # for k = 1:size(D1.V,1), i = 1:round(Int, l/2)
        #     D1.V[k,i], D1.V[k,l+1-i] = D1.V[k,l+1-i], D1.V[k,i]
        # end
        # arrange the values of S in non-decreasing order
        s = reverse(D1.S)
        # swap S
        S = zeros(Float64, p, l)
        j = 1
        for i = l-q2+1:l
            S[j, i] = s[j]
            j += 1
        end

        # set alpha(1:l-q2) = 1, beta(l-q2) = 0
        alpha = fill(0.0, l)
        beta = fill(0.0, l)
        if l-q2 >= 1
            for i  = 1:l-q2
                alpha[i] = 1
                beta[i] = 0
            end
        end

        j = 1
        for i = l-q2+1:l
            beta[i] = s[j]
            j += 1
        end

        # find r such that beta(r) <= 1/sqrt(2) <=beta(r+1)
        r = 1
        for i = l-q2+1:l-1
            if beta[i] <= sqrt(0.5) && sqrt(0.5) <= beta[i+1]
                r = i
                break
            end
        end

        # compute T = Q1Z
        T = Q1 * D1.V
        row, col = size(T)

        # compute QR decomposition of T
        D2 = qr(T)

        Q_ = D2.Q * Matrix(I, m, m)

        # sanitize R so that it only has non-negative diagonal entries
        dia = sign.(diag(D2.R))
        R = Diagonal(dia) * D2.R
        # In the following case, R is compactly stored as l by l matrix,
		# we need to resize it to m by l
        if m > l
            s = fill(1.0, m)
            for i = 1:length(dia)
                s[i] = dia[i]
            end
            R = [R; zeros(Float64, m-l, l)]
            U = Q_ * Diagonal(s)
        else
            U = Q_ * Diagonal(dia)
        end
        # println("Time for sanitization: " + time1 + "s.")

        #R2 = R[l-q2+1:r, l-q2+1:r]
        @views R2 = R[l-q2+1:r, l-q2+1:r]
        @views R3 = R[r+1:q1, r+1:l]

        # compute the SVD of R3
        D3 = svd(R3, full = true)

        # set alpha(r+1:q1)
        j = 1
        for i = r+1:q1
            alpha[i] = D3.S[j]
            j += 1
        end

        # set alpha(q1+1:l) = 0, beta(q1+1:l)=1
        if q1+1 <=l
            for i = q1+1:l
                alpha[i] = 0
                beta[i] = 1
            end
        end

        # set alpha(l-q2+1:r) = diag(R2)
        dia = diag(R2)
        j = 1
        for i = l-q2+1:r
            alpha[i] = dia[j]
            j = j + 1
        end

        # update U
        @views U[:,r+1:q1] = U[:,r+1:q1]*D3.U

        # update Z
        @views D1.V[:,r+1:l] = D1.V[:,r+1:l] * D3.V

        # set W
        @views S_ = Matrix(Diagonal(beta[r+1:q2]))
        @views W = S_ * D3.V[1:q2-r,1:q2-r]

        # compute QR of W
        D4 = qr(W)
        row = size(W)[1]
        col = size(W)[2]
        Rq = D4.Q * Matrix(I, row, row)

        # sanitize Rw so that it only has non-negative diagonal entries
        dia = sign.(diag(D4.R))
        Rw = Diagonal(sign.(diag(D4.R))) * D4.R
        # In the following case, Rw is compactly stored as col by col matrix,
		# we need to resize it to row by col
        if  row > col
            s = fill(1.0, row)
            for i = 1:length(dia)
                s[i] = dia[i]
            end
            Rw = [Rw; zeros(Float64, row-col, col)]
            # Rw = vcat(Rw, zeros(Float64, row-col, col))
            Uw = Rq * Diagonal(s)
        else
            Uw = Rq * Diagonal(dia)
        end
        # println("Time for sanitization: " + time2 + "s.")


        # update V
        n = min(r, l-q2)
        @views D1.U[:,1+r-n:q2-n] = D1.U[:,1+r-n:q2-n] * Uw

        # pack return value
        # V = D1.U
        # Z = D1.V
        C = zeros(Float64, m, l)
        if m < l
            for i = 1:m
                @views C[i,i] = alpha[i]
            end
        else
            for i = 1:l
                @views C[i,i] = alpha[i]
            end
        end

        if option == 1
            return U, D1.U, D1.V, C, S
        else
            return U, D1.U, D1.V, alpha, beta
        end

    # case 2: m > p
    else
        # compute the full SVD of Q1
        D1 = svd(Q1, full = true)

        # set alpha(q1+1:l) = 0, beta(q1+1:l) = 0
        alpha = fill(0.0, l)
        beta = fill(0.0, l)
        if l >= q1+1
            for i = q1+1:l
                alpha[i] = 0
                beta[i] = 1
            end
        end

         for i = 1:q1
            alpha[i] = D1.S[i]
        end

        # find r such that 1>=alpha(1)>=...>=alpha(r)>=1/sqrt(2)
		# >alpha(r+1)>=...>=alpha(q1)>=0
        k_p = max(l-p, 0)
        r = 0
        if alpha[q1] > sqrt(0.5)
            r = q1 - k_p
        else
            j = k_p + 1
            for i = j:q1
                if alpha[i] < sqrt(0.5)
                    r = i - j
                    break
                end
            end
        end

        # compute T = Q2 * Z
        Z = D1.V
        T = Q2 * Z

        # compute the QL decomposition of T
        A, tau = LAPACK.geqlf!(T)
        A_ = copy(A)
        Q = LAPACK.orgql!(A, tau, length(tau))
        # @time Q = QL_generateQ(A, tau)

        # if p >= l
        #     @views L = LowerTriangular(A[p-l+1:p,1:l])
        #     # sanitize L so that it only has non-negative diagonal entries
        #     dia = Diagonal(sign.(diag(L)))
        #     L = dia* L
        #     @views Q[1:p, p-l+1:p] = Q[1:p, p-l+1:p] * dia
        #     L_= [zeros(Float64, p-l, l); L]
        # else
        #     @views L = LowerTriangular(A[1:p,l-p+1:l])
        #     # sanitize L so that it only has non-negative diagonal entries
        #     dia = Diagonal(sign.(diag(L)))
        #     L = dia* L
        #     @views A_ = dia*A[1:p,1:l-p]
        #     Q = Q*dia
        #     L_= [A_ L]
        # end

        if p >= l
            @views L = LowerTriangular(A_[p-l+1:p,1:l])
            # sanitize L so that it only has non-negative diagonal entries
            dia = Diagonal(sign.(diag(L)))
            L = dia* L
            @views Q[1:p, p-l+1:p] = Q[1:p, p-l+1:p] * dia
            L_= [zeros(Float64, p-l, l); L]
        else
            @views L = LowerTriangular(A_[1:p,l-p+1:l])
            # sanitize L so that it only has non-negative diagonal entries
            dia = Diagonal(sign.(diag(L)))
            L = dia* L
            @views A_ = dia*A_[1:p,1:l-p]
            Q = Q*dia
            L_= [A_ L]
        end

        @views L1 = L_[p-q2+1:p-q2+r,1:l-q2+r]
        @views L2 = L_[p-q2+r+1:q1-l+p,r+1+l-q2:q1]
        # L1 = L_[p-q2+1:p-q2+r,1:l-q2+r]
        # L2 = L_[p-q2+r+1:q1-l+p,r+1+l-q2:q1]

        # compute the SVD of L1
        D2 = svd(L1, full = true)
        Vl = D2.U
        Zl = D2.V

        # arrange the values of S in non-decreasing order
        s = reverse(D2.S)
        # swap beta
        j = 1
        for i = l-q2+1:r+l-q2
            beta[i] = s[j]
            j += 1
        end

        # swap Vl
        @views Vl[:,1:r] = Vl[:,r:-1:1]

        # swap Zl
        @views Zl[:,1:r+l-q2] = Zl[:,r+l-q2:-1:1]

        # set alpha(1:l-q2) = 1, beta (1:l-q2) = 0
        if l-q2 >= 1
            for i = 1: l-q2
                alpha[i] = 1
                beta[i] = 0
            end
        end

        # set beta(r+l-q2+1:q1) = diag(L2)
        dia = diag(L2)
        j = 1
        for i = r+l-q2+1:q1
            beta[i] = dia[j]
            j += 1
        end

        # set V_
        V_ = zeros(Float64, p, p)
        V_[p-q2+1:p-q2+r,1:r] = Vl
        V_[1:p-q2,q2+1:p] =  Matrix{Float64}(I, p-q2, p-q2)
        V_[p-q2+r+1:p,r+1:q2] = Matrix{Float64}(I, q2-r, q2-r)

        # update V
        V = Q*V_

        # update Z
        @views Z[1:l, 1:r+l-q2] = Z[1:l, 1:r+l-q2] * Zl

        # set W
        @views S_ = Matrix(Diagonal(alpha[1:r+l-q2]))
        W = S_ * Zl

        # compute QR of W
        D4 = qr(W)
        row = size(W)[1]
        col = size(W)[2]
        Qw_ = D4.Q * Matrix(I, row, row)
        # sanitize Rw so that it only has non-negative diagonal entries
        Rw = Diagonal(sign.(diag(D4.R))) * D4.R
        # In the following case, Rw is compactly stored as col-by-col matrix,
		# we need to resize it to row-by-col
        dia = sign.(diag(D4.R))
        if  row > col
            s = fill(1.0, row)
            for i = 1:length(dia)
                s[i] = dia[i]
            end
            Rw = [Rw; zeros(Float64, row-col, col)]
            Qw = Qw_ * Diagonal(s)
        else
            Qw = Qw_ * Diagonal(dia)
        end

        # update U
        @views D1.U[1:m,1:r+l-q2] = D1.U[1:m,1:r+l-q2] * Qw

        # formulate C
        C = zeros(Float64, m, l)
        for i = 1:q1
            C[i, i] = alpha[i]
        end

        # formulate S
        S = zeros(Float64, p, l)
        if p >= l
            for i = 1:l
                S[i, i] = beta[i]
            end
        else
            j = 1
            for i = l-p+1:l
                S[j, i] = beta[i]
                j += 1
            end
        end
        if option == 1
            return D1.U, V, Z, C, S
        else
            return D1.U, V, Z, alpha, beta
        end
    end
end
