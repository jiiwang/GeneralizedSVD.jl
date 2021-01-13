"""
and stability of CS Decomposition

First, we generate a m+p by m+p random matrix A
Then, we compute QR decomposition of A, A = QR
Afterwards, we choose Q1 = Q[1:m, 1:l]
Q2 = Q[m+1:m+p, 1:l], such that:

      Q1' * Q1 + Q2' * Q2 = I

Call csd function, return U, V, Z, C, S.
Calculate and print time elapsed for computing CSD.

Assume ||Q1' * Q1 + Q2' * Q2 - I|| ≈ eps,
compute:

      res_q1 = ||U * C * Z' - Q1||/(max(m,l)* ||Q1|| * eps)

      res_q2 = ||V * S * Z' - Q2||/(max(p,l)* ||Q2|| * eps)

      orthog_cs = ||C' * C + S' * S -I||/(max(m,p,l) * eps)

      orthog_u = ||U' * U - I||/(m * eps)

      orthog_v = ||V' * V - I||/(p * eps)

      orthog_z = ||Z' * Z - I||/(l * eps)

argument
m: num of rows in Q1
p: num of rows in Q2
l: num of cols in Q1, Q2

no return
"""

function test(m, p, l)
    if m + p < l
        println("input error, m + p >= l is required.")
        return
    else
        F = qr(randn(m+p, l))
        # Q1 = view(Q_, 1:m, 1:l)
        # Q2 = view(Q_, m+1:m+p, 1:l)

        Q1 = Matrix(F.Q)[1:m, :]
        Q2 = Matrix(F.Q)[m+1:m+p, :]

        e = eps(Float64)

        @time U, V, Z, C, S = csd(Q1, Q2, 1)

        res_q1 = opnorm(U * C * Z' - Q1, 1)/(max(m,l)*opnorm(Q1, 1)*e)
        res_q2 = opnorm(V * S * Z' - Q2, 1)/(max(p,l)*opnorm(Q2, 1)*e)
        orthog_cs = opnorm(C' * C + S' * S - I, 1)/(max(m,p,l)*e)
        orthog_u = opnorm(U' * U - I, 1)/(m*e)
        orthog_v = opnorm(V' * V - I, 1)/(p*e)
        orthog_z = opnorm(Z' * Z - I, 1)/(l*e)
        
        println("res_q1: ", res_q1)
        println("res_q2: ", res_q2)
        println("orthog_cs: ", orthog_cs)
        println("orthog_u: ", orthog_u)
        println("orthog_v: ", orthog_v)
        println("orthog_z: ", orthog_z)
    end
end

# This function works similarly to the previous test function,
# Instead of evaluating time performance and stability,
# it returns the computed C, S.
#
# argument
# m: num of rows in Q1
# p: num of rows in Q2
# l: num of cols in Q1/Q2
#
# returns C, S
function testCS(m, p, l)
    if m + p < l
        println("input error, m + p >=l is required.")
        return
    else
        Q, R = qr(randn(m+p, m+p))
        Q_ = Q*Matrix(I,m+p,m+p)
        Q1 = view(Q_, 1:m, 1:l)
        Q2 = view(Q_, m+1:m+p, 1:l)
        U, V, Z, C, S = csd(Q1, Q2, 1)
        return C, S
    end
end

# This function works similarly to the previous test function,
# Instead of evaluating time performance and stability,
# it returns the computed U, V, Z, C, S.
#
# argument
# m: num of rows in Q1
# p: num of rows in Q2
# l: num of cols in Q1/Q2
#
# returns U, V, Z, C, S
function testProduct(m, p, l)
    if m + p < l
        println("input error, m + p >=l is required.")
        return
    else
        Q, R = qr(randn(m+p, m+p))
        Q_ = Q*Matrix(I,m+p,m+p)
        Q1 = view(Q_, 1:m, 1:l)
        Q2 = view(Q_, m+1:m+p, 1:l)
        return U, V, Z, C, S = csd(Q1, Q2, 1)
    end
end

# Example 1 in Chapter 3 of Jenny's thesis
function test1a()
    Q1 = [1/sqrt(7) 0 0 1/sqrt(3);
          1/sqrt(7) -2/sqrt(10) -1/2 -1/(2*sqrt(3));
          1/sqrt(7) -1/sqrt(10) 1/4  3/(4*sqrt(3));
          1/sqrt(7) 0 3/4 -3/(4*sqrt(3));
          1/sqrt(7) 0 0 0]
    Q2 = [1/sqrt(7) 1/sqrt(10) -1/4 -3/(4*sqrt(3));
          1/sqrt(7) 2/sqrt(10) -1/4 1/(4*sqrt(3))]
    U, V, Z, alpha, beta = csd(Q1, Q2, 0)
end

function test1b()
    Q1 = [1/sqrt(7) 0 0 1/sqrt(3);
          1/sqrt(7) -2/sqrt(10) -1/2 -1/(2*sqrt(3));
          1/sqrt(7) -1/sqrt(10) 1/4  3/(4*sqrt(3));
          1/sqrt(7) 0 3/4 -3/(4*sqrt(3));
          1/sqrt(7) 0 0 0]
    Q2 = [1/sqrt(7) 1/sqrt(10) -1/4 -3/(4*sqrt(3));
          1/sqrt(7) 2/sqrt(10) -1/4 1/(4*sqrt(3))]
    U, V, Z, C, S = csd(Q1, Q2, 1)
end

# Example 2 in Chapter 3 of Jenny's thesis
function test2()
    Q1 = [1/sqrt(7) 0 0 1/sqrt(3);
          1/sqrt(7) -2/sqrt(10) -1/2 -1/(2*sqrt(3));
          1/sqrt(7) -1/sqrt(10) 1/4  3/(4*sqrt(3))]
    Q2 = [1/sqrt(7) 0 3/4 -3/(4*sqrt(3));
          1/sqrt(7) 0 0 0;
          1/sqrt(7) 1/sqrt(10) -1/4 -3/(4*sqrt(3));
          1/sqrt(7) 2/sqrt(10) -1/4 1/(4*sqrt(3))]
    U, V, Z, C, S = csd(Q1, Q2, 1)
end
