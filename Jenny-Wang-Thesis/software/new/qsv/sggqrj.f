      SUBROUTINE SGGQRJ  ( M, N, A, LDA, B, LDB, U, LDU, V, LDV,
     $                     WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 0.1) --
*     Univ. of California, Davis
*     September, 2005
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDU, LDV, M, N, LWORK, INFO
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ),
     $                   U( LDU, * ), V( LDV, * ),
     $                   WORK( LWORK )
*     ..
*
*
* PURPOSE:
* ========
*
* SGGQRJ computes the compact QR factorization of a (M+N)-by-N
* matrix where A is M-by-N upper-trapezoidal and B is
* N-by-N upper-triangular:
*
*      Q R = (A' B')'.
*
*
* Arguments
* =========
*
* M       (input) INTEGER
*         The number of rows of the matrix A.  M >= 0.
*
* N       (input) INTEGER
*         The number of columns of the matrix A
*         and he number of rows and columns of the matrix B.
*         N >= 0, N >= M.
*
* A       (input/output) REAL array, dimension (LDA, N)
*         On entry, the M-by-N upper-trapezoidal matrix A.
*         On exit, the first min( M, N ) rows of the upper triangular
*         matrix R.
*
* LDA     (input) INTEGER
*         The leading dimension of the array A. LDA >= max(1, M).
*
* B       (input/output) REAL array, dimension (LDB, N)
*         On entry, the N-by-N upper triangular matrix B
*         On exit, the origional content of B is destroyed.
*         If M < N, B contains the M+1 to N rows of the upper
*         triangular matrix R.
*
* LDB     (input) INTEGER
*         The leading dimension of the array B. LDB >= max(1, N).
*
* U       (output) REAL array, dimension (LDU, N)
*         U contains the first M rows of the orthogonal matrix Q.
*
* LDU     (input) INTEGER
*         The leading dimension of the array U. LDU >= max(1, M).
*
* V       (output) REAL array, dimension (LDV, N)
*         V contains the last N rows of the orthogonal matrix Q.
*
* LDV     (input) INTEGER
*         The leading dimension of the array V. LDV >= max(1, N).
*
* WORK    (workspace) REAL array, dimension (LWORK)
*
* LWORK   (input) INTEGEER
*         The dimension of the array WORK.
*	  LWORK = M+N.
*
* INFO    (output) INTEGER
*         = 0:  successful exit
*         < 0:  if INFO = -i, the i-th argument had an illegal value.
*
* =====================================================================
*
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K
      REAL               C, S, R
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASET, SLARTG, SROT, SCOPY
*     ..
*     .. Intrinsic Functions
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      I = MAX( 1, M )
      J = Max( 1, N )
      INFO = 0
      IF( ( M.LT.0 ) .OR. ( M.GT.N ) ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.I ) THEN
        INFO = -4
      ELSE IF( LDB.LT.J ) THEN
        INFO = -6
      ELSE IF( LDU.LT.I ) THEN
        INFO = -8
      ELSE IF( LDV.LT.J ) THEN
        INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'SGGQRJ', -INFO )
        RETURN
      END IF
*
*     Executable statements
*
*--------
*     STEP 1. Initialize U and V so that (U' V')' is the identity
*--------
*
      CALL SLASET( 'Full', M, N, ZERO, ONE, U, LDU  )
      CALL SLASET( 'Full', N, N, ZERO, ZERO, V, LDV )
      DO I = 1, N - M
         V( I, I + M ) = ONE
      END DO
*--------
*     STEP 2.
*--------
      DO I = 1, N
         CALL SLASET( 'Full', M+N, 1, ZERO, ZERO, WORK, 1 )
         WORK( M + I ) = ONE
         K = MIN( N, M - 1 + I )
         DO J = I, K
*--------
*     STEP 2.1.  Eliminate elements in B
*--------
            IF( J.LE.M ) THEN
*             use the j-th row of A to eliminate the i-th row of B
*
              CALL SLARTG( A( J, J ), B( I, J ), C, S, R )
              A( J, J ) = R
              B( I, J ) = ZERO
              IF( J+1 .LE. N ) THEN
*               the rest of the J+1:N columns are effected accordingly
*
                CALL SROT( N-J, A( J,J+1 ), LDA, B( I, J+1 ), LDB,
     $                     C, S)
              END IF
*
            ELSE
*             use the (j-m)-th row of B to eliminate the i-th row of B
*
              CALL SLARTG( B( J-M, J ), B( I, J ), C, S, R )
              B( J-M, J ) = R
              B(   I, J ) = ZERO
*
              IF( J+1 .LE. N ) THEN
*               the rets of the J+1:N columns are affected accordingly
*
                CALL SROT( N-J, B( J-M, J+1 ), LDB, B( I, J+1 ), LDB,
     $                     C, S)
              END IF
            END IF
*
*--------
*     STEP 2.2.  Update two columns in A and B
*--------
*
            CALL SROT( M, U( 1, J ), 1, WORK( 1 ), 1,  C, S )
            CALL SROT( N, V( 1, J ), 1, WORK( M+1 ), 1, C, S )
         END DO
*
*--------
*     STEP 2.3.  Place contents of WORK(M+N) to U and V if needed
*--------
*
         IF( M+I .LE. N ) THEN
           CALL SCOPY( M, WORK       , 1, U(1, M+I), 1 )
           CALL SCOPY( N, WORK( M+1 ), 1, V(1, M+I), 1 )
         END IF
      END DO
*
      RETURN
*
*     End of SGGQRJ
*
      END
