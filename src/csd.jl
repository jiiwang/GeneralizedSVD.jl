"""
	csd2by1(Q1, Q2) computes the the Cosine-sine decomposition of
	an (m+l) by l orthogonal matrix Q = [Q1; Q2] where m <= l, such that:

```math
	Q1 = U * C * Z',    Q2 = V * S * Z'
```

- U is m-by-m orthogonal matrix,
- V is l-by-l orthogonal matrix,
- Z is l-by-l orthogonal matrix,
- C and S are  m-by-l and l-by-l ``diagonal" matrices such that
C' * C + S' * S = I. C and S have the following structures
(both alpha and beta are of length l):

(1) If m == l:
        C =  diag(alpha(1), ... , alpha(l))

        S =  diag(beta(1), ... , beta(l))

(2) If m < l:
        C =                       m                 l-m
              m ( diag(alpha(1), ... , alpha(m))     0 )

        S =                       m                 l-m
              m ( diag(beta(1), ... , beta(m)))      0 )
            l-m (                 0                  I )

# arguments:
- Q1: m by l matrix,
- Q2: l by l matrix.

# returns:
- U is m-by-m matrix,
- V is l-by-l matrix,
- Z is l-by-l matrix,
- cosine is an array of length l,
- sine is an array of length l.
"""
# csd2by1() calls csd2by1!() before making a copy of input matrices A and B,
# thus A, B will not be overwritten
function csd2by1(Q1::StridedMatrix{T}, Q2::StridedMatrix{T}) where T
	return csd2by1!(copy(Q1),copy(Q2))
end

function csd2by1!(Q1::StridedMatrix{T}, Q2::StridedMatrix{T}) where T
    m, l = size(Q1)
	if size(Q2)[1] != size(Q2)[2]
		throw(ArgumentError("Q2's row number must equal its column number"))
	elseif m > size(Q2)[1]
		throw(ArgumentError("Q1's row number must be less than that of Q2"))
	elseif l != size(Q2)[2]
		throw(ArgumentError("Q1's column number must equal that of Q2"))
	end

    # Step 1: compute the SVD of Q2
	# May cause LAPACK Exception: S/DBDSDC did not converge
	D1 = svd!(Q2)

	# Step 2.1: rearrange the values of S in non-decreasing order
	sine = reverse(D1.S)

    # Step 2.2: reverse V's col
	#     Assign array slicing to D1.U will throw an error:
	#     D1.U = @views D1.U[:, l:-1:1]
	#     ERROR: setfield! immutable struct of type SVD cannot be changed
	D1.U[:,1:l] = @views D1.U[:, l:-1:1]

    # Step 2.3: reverse Z's col
    D1.V[:,1:l] = @views D1.V[:,l:-1:1]

    # Step 3: find r such that sine(r) <= 1/sqrt(2) < sine(r+1),
	# binary search is used
	r = binSearch(sine)

    # Step 4: compute X = Q1Z
    X = Q1 * D1.V

    # Step 5: compute QR decomposition of X
    D2 = qr!(X)

    # Step 5a: sanitize R so that it only has non-negative diagonal entries
	R = D2.R
    dia1 = sign.(diag(R))
	# update dia1 if R has zeros in the diagonal
	dia1[dia1 .== T(0.0)] .= T(1.0)
	cosine = zeros(T, l)
	@views cosine[1:r] = abs.(diag(R)[1:r])
    R = Diagonal(dia1) * R
	U = Matrix{T}(D2.Q) * Diagonal(dia1)

	# Special case: exit
	if m <= r
		return U, D1.U, D1.V, cosine, sine
	end

	# Step 6: refine R
	R22 = @views R[r+1:m,r+1:l]

	# Step 6a: compute the full SVD of R22
	D3 = svd!(R22, full = true)

	# Step 6b: update the (r+1)-th to m-th cols of U
	U[:,r+1:m] = @views U[:,r+1:m]*D3.U

	# Step 6c: update the last (l-r) cols of Z
	@views D1.V[:,r+1:l] = D1.V[:,r+1:l] * D3.V

	# Step 6d: update cosine and sine
	@views sine[m+1:l] = ones(T, l-m)
	@views cosine[r+1:m] = D3.S

	# Step 7: refine V
    # Step 7a: set W
	W = @views Diagonal(sine[r+1:l])* D3.V

    # Step 7b: compute QR of W
    D4 = qr!(W)
    Qw = D4.Q * Matrix{T}(I, l-r, l-r)

    # Step 7c: sanitize D4.R so that it only has non-negative diagonal entries
    dia2 = sign.(diag(D4.R))
    Qw = Qw * Diagonal(dia2)

    # Step 7d: update V
    @views D1.U[:,r+1:l] = D1.U[:,r+1:l] * Qw

    return U, D1.U, D1.V, cosine, sine
end

# binSearch() modifies standard binary search
# to find r such that beta(r) <= 1/sqrt(2) < beta(r+1),
# return 1 if no such r exist.
function binSearch(a::Vector{T}) where T
	r = 1
	low, high = 1, length(a)
	thld = sqrt(T(0.5))
	while low <= high
		mid = floor(Integer, (high - low)/2) + low
		r = mid
		if mid+1 <= length(a) && a[mid+1] <= thld
			low = mid + 1
		elseif a[mid] > thld
			high = mid - 1
		else
			return mid
		end
	end
	return r
end
