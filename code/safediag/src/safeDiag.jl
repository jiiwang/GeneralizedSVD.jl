# module CSD

using LinearAlgebra

# export csd

"""
	safeDiag(Q1, Q2)

Compute the safe diagonalization (a special case of the CSD) of an
(m+l) by l orthogonal matrix Q = (Q1' Q2')' where m <= l, such that:

Q1 = U  C  Z',  Q2 = V  S  Z'

and return U, V, Z, C, S in order.

Specifically, U, V and Z are orthogonal matrices of size m-by-m, l-by-l,
and l-by-l, respectively.
C and S are  m-by-l and l-by-l "diagonal" matrices such that
C' * C + S' * S = I.

C and S have the following structures (both alpha and beta are of length l):

(1) If m == l:
        C =  diag(alpha(1), ... , alpha(l))

        S =  diag(beta(1), ... , beta(l))

(2) If m < l:
        C =                       m                 l-m
              m ( diag(alpha(1), ... , alpha(m))     0 )

        S =                       m                 l-m
              m ( diag(beta(1), ... , beta(m)))      0 )
            l-m (                 0                  I )


# Arguments
- Q1: m by l matrix
- Q2: l by l matrix

"""

function safeDiag(Q1, Q2)
    m, l = size(Q1)
	if size(Q2)[1] != size(Q2)[2]
		throw(ArgumentError("Q2's row number must equal its column number"))
	elseif m > size(Q2)[1]
		throw(ArgumentError("Q1's row number must be less than that of Q2"))
	elseif l != size(Q2)[2]
		throw(ArgumentError("Q1's column number must equal that of Q2"))
	end

    # Step 1: compute the SVD of Q2
    D1 = svd(Q2)

	# Step 2.1: rearrange the values of S in non-decreasing order
	beta = reverse(D1.S)
    S = Diagonal(beta)

    # Step 2.2: reverse V's col
	# Comment:
	#     Assign array slicing to D1.U will throw an error:
	#     D1.U = @views D1.U[:, l:-1:1]
	#     ERROR: setfield! immutable struct of type SVD cannot be changed
	D1.U[:,1:l] = @views D1.U[:, l:-1:1]

    # Step 2.3: reverse Z's col
    D1.V[:,1:l] = @views D1.V[:,l:-1:1]

	# check residual of Q2
	# println("res_q2 after Step 2: ", opnorm(D1.U * S * (D1.V)' - Q2, 1)/(l*opnorm(Q2, 1)*eps(Float64)))

    # Step 3: find r such that beta(r) <= 1/sqrt(2) < beta(r+1),
	# binary search is used
	r = binSearch(beta)

    # Step 4: compute T = Q1Z
    T = Q1 * D1.V

    # Step 5: compute QR decomposition of T
    D2 = qr(T)
    U = Matrix(D2.Q)

    # Step 5a: sanitize R so that it only has non-negative diagonal entries
    dia1 = sign.(diag(D2.R))
    R = Diagonal(dia1) * D2.R

	U = U * Diagonal(dia1)

	# Check residual of Q1
	# println("res_q1 after Step 5: ", opnorm(U * R * (D1.V)' - Q1, 1)/(l*opnorm(Q1, 1)*eps(Float64)))

	# Special case: exit
	if (m == l && l <= r) || (m < l && m <= r)
		return U, D1.U. D1.V, R, S
	end

	# Step 6: refine R
	# Two cases can be combined as one procedure
	R22 = @views R[r+1:m,r+1:l]

	# Step 6a: compute the SVD of R22
	D3 = svd(R22, full = true)

	# Step 6b: update the (r+1)-th to m-th cols of U
	@views U[:,r+1:m] = U[:,r+1:m]*D3.U

	# Step 6c: update the last (l-r) cols of Z
	@views D1.V[:,r+1:l] = D1.V[:,r+1:l] * D3.V

	# Step 6d: pack R to formulate C
	C = zeros(Float64, m, l)
    for i = 1:r
        @views C[i,i] = R[i,i]
    end
	for i = 1:m-r
		@views C[i+r,i+r] = D3.S[i]
	end

	# Check residual of Q1
	# println("res_q1 after Step 6: ", opnorm(U * C * (D1.V)' - Q1, 1)/(l*opnorm(Q1, 1)*eps(Float64)))

	# Step 7: refine V
    # Step 7a: set W
    W = @views S[r+1:l,r+1:l] * D3.V

    # Step 7b: compute QR of W
    D4 = qr(W)
    Qw = D4.Q * Matrix(I, l-r, l-r)

    # Step 7c: sanitize D4.R so that it only has non-negative diagonal entries
    dia2 = sign.(diag(D4.R))
    Qw = Qw * Diagonal(dia2)

    # Step 7d: update V
    @views D1.U[:,r+1:l] = D1.U[:,r+1:l] * Qw

	# Check residual of Q1
	# println("res_q2 after Step 7: ", opnorm(D1.U * S * (D1.V)' - Q2, 1)/(l*opnorm(Q2, 1)*eps(Float64)))

    return U, D1.U, D1.V, C, S
end

function binSearch(beta)
	r = 1
	low, high = 1, length(beta)
	thld = sqrt(0.5)
	while low <= high
		mid = floor(Integer, (high - low)/2) + low
		r = mid
		if beta[mid+1] <= thld
			low = mid + 1
		elseif beta[mid] > thld
			high = mid - 1
		else
			return mid
		end
	end
	return r
end
