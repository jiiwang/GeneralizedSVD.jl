"""
    householderqr(R1, R2) computes the QR decompsition of
    (m+l)-by-l block matrix [R1;R2]

```math
          l               l
    m  (  R1  )  =  m  (  Q1  ) * R
    l  (  R2  )     l  (  Q2  )
```

- Q1 is m-by-l orthogonal matrix,
- Q2 is l-by-l orthogonal matrix,
- R is l-by-l upper triangular matrix.

# arguments:
R1: m-by-l Float 32 or Float64 matrix,
R2: m-by-l Float 32 or Float64 matrix.

# returns:
Q1: m-by-l matrix,
Q2: l-by-l matrix,
R: l-by-l matrix.
"""
function householderqr(R1, R2)
    r, c = size(R1)
    F = qr!([R1;R2])
    Q = Matrix(F.Q)
    Q1 = Q[1:r,:]
    Q2 = Q[r+1:end,:]
    return Q1, Q2, F.R
end
