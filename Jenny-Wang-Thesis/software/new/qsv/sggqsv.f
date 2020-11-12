      SUBROUTINE SGGQSV( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B,
     $                   LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ,
     $                   IWORK, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 0.1) --
*     Univ. of California Davis,
*     September, 2005
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          JOBQ, JOBU, JOBV
      INTEGER            M, N, P, K, L, LDA, LDB, LDU, LDV, LDQ,
     $                   LWORK, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), Q( LDQ, * ), U( LDU, * ),
     $                   V( LDV, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SGGQSV computes the quotient ( or generalized ) singular value
*  decomposition (QSVD) of an M-by-N real matrix A and P-by-N
*  real matrix B:
*
*      U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R )
*
*  where U, V and Q are orthogonal matrices, and Z' is the transpose
*  of Z.  Let K+L = the effective numerical rank of the matrix (A',B')',
*  then R is a (K+L)-by-(K+L) nonsingular upper triangular matrix,
*  D1 and D2 are M-by-(K+L) and P-by-(K+L) "diagonal" matrices and
*  have one of the following structures:
*
*  If M >= K+L,
*
*         D1 =      K   L
*             K  (  I   0  )
*             L  (  0   C  )
*         M-K-L  (  0   0  )
*
*         D2 =      K   L
*             L  (  0   S  )
*           P-L  (  0   0  )
*
*                  N-K-L  K    L
*    ( 0 R ) = K (  0    R11  R12 )
*              L (  0     0   R22 )
*
*  where
*
*       C = diag( ALPHA( K+1 ), ... , ALPHA( K+L ) ),
*       S = diag( BETA( K+1 ), ..., BETA( K+L ) ),
*       C**2 + S**2 = I
*
*       R is stored in A( 1:K+L , N-K-L+1:N ) on exit.
*
*  (2) If M < K+L,
*
*        D1 =         K   M-K  K+L-M
*                 K ( I    0     0  )
*               M-K ( 0    C     0  )
*
*         D2 = 	      K  M-K   K+L-M
*               M-K ( 0   S     0  )
*             K+L-M ( 0   0     I  )
*             P-K-L ( 0   0     0  )
*
*                   N-K-L   K   M-K  K+L-M
*    ( 0 R ) =    K ( 0    R11  R12  R13  )
*               M-K ( 0     0   R22  R23  )
*             K+L-M ( 0     0    0   R33  )
*
*  where
*
*       C = diag( ALPHA(K+1), ... , ALPHA(M) ),
* 	    S = diag( BETA(K+1), ..., BETA(M) ),
*       C**2 + S**2 = I
*
*    (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored
*    ( 0  R22 R23 )
*    in B(1:K+L-M, N+M-K-L+1:N) on exit.
*
*  The routine computes C, S, R, and optionally the orthogonal
*  transformation matrices U, V and Q.
*
*  In particular, if B is an N-by-N nonsingular matrix, then the GSVD of
*  A and B implicitly gives the SVD of A*inv(B):
*                       A*inv(B) = U*(D1*inv(D2))*V'.
*  If ( A',B')' has orthonormal columns, then the GSVD of A and B is
*  also equal to the CS decomposition of A and B. Furthermore, the GSVD
*  can be used to derive the solution of the eigenvalue problem:
*                       A'*A x = lambda* B'*B x.
*  In some literature, the GSVD of A and B is presented in the form
*                   U'*A*X = ( 0 D1 ),   V'*B*X = ( 0 D2 )
*  where U and V are orthogonal and X is nonsingular, D1 and D2 are
*  ``diagonal''.  The former GSVD form can be converted to the latter
*  form by taking the nonsingular matrix X as
*
*                       X = Q*( I   0    )
*                             ( 0 inv(R) ).
*
*  Arguments
*  =========
*
*  JOBU    (input) CHARACTER*1
*          = 'U':  Orthogonal matrix U is computed;
*          = 'N':  U is not computed
*
*  JOBV    (input) CHARACTER*1
*          = 'V':  Orthogonal matrix V is computed;
*          = 'N':  V is not computed
*
*  JOBQ    (input) CHARACTER*1
*          = 'Q':  Orthogonal matrix Q is computed;
*          = 'N':  Q is not computed
*
*  NOTE: It is assumed that JOBU, JOBV, JOBQ are all turned on
*        or off at the same time.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B.  N >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
*  K       (output) INTEGER
*  L       (output) INTEGER
*          On exit, K and L specify the dimension of the subblocks
*          described in the Purpose section.
*          K + L = effective numerical rank of (A',B')'.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A contains the triangular matrix R, or part of R.
*          See Purpose for details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input/output) REAL array, dimension (LDB,N)
*          On entry, the P-by-N matrix B.
*          On exit, B contains the triangular matrix R if M-K-L < 0.
*          See Purpose for details.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDA >= max(1,P).
*
*  ALPHA   (output) REAL array, dimension (N)
*  BETA    (output) REAL array, dimension (N)
*          On exit, ALPHA and BETA contain the generalized singular
*          value pairs of A and B;
*            ALPHA( 1:K ) = 1,
*            BETA( 1:K )  = 0,
*          and if M-K-L >= 0,
*            ALPHA( K+1:K+L ) = C,
*            BETA( K+1:K+L )  = S,
*          or if M-K-L < 0,
*            ALPHA( K+1:M )=C, ALPHA( M+1:K+L )=0
*            BETA( K+1:M ) =S, BETA( M+1:K+L ) =1
*          and
*            ALPHA(K+L+1:N) = 0
*            BETA(K+L+1:N)  = 0
*
*  U       (output) REAL array, dimension (LDU,M)
*          If JOBU = 'U', U contains the M-by-M orthogonal matrix U.
*          If JOBU = 'N', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U. LDU >= max(1,M) if
*          JOBU = 'U'; LDU >= 1 otherwise.
*
*  V       (output) REAL array, dimension (LDV,P)
*          If JOBV = 'V', V contains the P-by-P orthogonal matrix V.
*          If JOBV = 'N', V is not referenced.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V. LDV >= max(1,P) if
*          JOBV = 'V'; LDV >= 1 otherwise.
*
*  Q       (output) REAL array, dimension (LDQ,N)
*          If JOBQ = 'Q', Q contains the N-by-N orthogonal matrix Q.
*          If JOBQ = 'N', Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= max(1,N) if
*          JOBQ = 'Q'; LDQ >= 1 otherwise.
*
*  IWORK   (workspace) INTEGER array, dimension ( N )
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*
*  LWORK   (input) INTEGEER
*         The dimension of the array WORK.
*         LWORK >= 6 * MIN( P, N )**2 + 6 * N + max( 3M, N, P )
*
*  INFO    (output)INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Internal Parameters
*  ===================
*
*  TOLA    REAL
*  TOLB    REAL
*          TOLA and TOLB are the thresholds to determine the effective
*          rank of (A',B')'. Generally, they are set to
*                   TOLA = MAX(M,N)*norm(A)*MACHEPS,
*                   TOLB = MAX(P,N)*norm(B)*MACHEPS.
*          The size of TOLA and TOLB may affect the size of backward
*          errors of the decomposition.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ...
*     .. Local Scalars ..
      LOGICAL            WANTQ, WANTU, WANTV, WANTALL
      INTEGER            I, J, G, T, IWK, LWK,
     $                   LDQ1, LDQ2, LDWU, LDWV, LDWZ, IQ1, IQ2,
     $                   IU, IV, IZ
      REAL               ANORM, BNORM, TOLA, TOLB, ULP, UNFL
      REAL               S_PRE, S_QRJ, S_CSD, S_POS, S_ALL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANGE, SECOND
      EXTERNAL           LSAME, SLAMCH, SLANGE, SECOND
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SORCSD, SGGQRJ, SGGSVP,
     $                   SLACPY, SLASET, XERBLA, SSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..

*
*     Test the input parameters
*
      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
      WANTQ = LSAME( JOBQ, 'Q' )
*
      IF( WANTU .OR. WANTV .OR. WANTQ ) THEN
        WANTALL = 'Y'
      ELSE
        WANTALL = 'N'
      END IF
*
      S_PRE = ZERO
      S_QRJ = ZERO
      S_CSD = ZERO
      S_POS = ZERO
      S_ALL = ZERO
*
      INFO = 0
      IF( .NOT.( WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
         INFO = -16
      ELSE IF( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) THEN
         INFO = -18
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGGQSV', -INFO )
         RETURN
      END IF
      G = MIN( M+P, N )
*
*     Compute the Frobenius norm of matrices A and B
*
      ANORM = SLANGE( '1', M, N, A, LDA, WORK )
      BNORM = SLANGE( '1', P, N, B, LDB, WORK )
*
*     Get machine precision and set up threshold for determining
*     the effective numerical rank of the matrices A and B.
*
      ULP = SLAMCH( 'Precision' )
      UNFL = SLAMCH( 'Safe Minimum' )
      TOLA = MAX( M, N )*MAX( ANORM, UNFL )*ULP
      TOLB = MAX( P, N )*MAX( BNORM, UNFL )*ULP

*      print *, "A = ..."
*      call ptmatrix(M, N, A, LDA)
*      print *, "B = ..."
*      call ptmatrix(P, N, B, LDB)
*
*     1. Preprocessing
*
      S_ALL = SECOND()
      S_PRE = SECOND()
*
      CALL SGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA,
     $             TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK,
     $             WORK, WORK( N + 1 ), INFO )
      S_PRE = SECOND() - S_PRE
*
      T = MIN( L, M-K )
      LDQ1 = M+1
      LDQ2 = P+1
      LDWU  = T+1
      LDWV  = L+1
      LDWZ  = L+1
      IQ1  = 1
      IQ2  = LDQ1 * L + IQ1
      IU   = LDQ2 * L + IQ2
      IV   = LDWU  * T + IU
      IZ   = LDWV  * L + IV
      IWK   = LDWZ  * L + IZ
      LWK   = LWORK - IWK
*      print *, "% K = , TOLA", K, TOLA
*      print *, "% L = , TOLB", L, TOLB
*      pause
*
*      ---- CHECK Result ----
*      print *, "Matrix = []; % U = "
*      CALL printMatrix(M, M, U, LDU )
*      print *, "U = Matrix; "
*      print *, "Matrix = []; % V = "
*      CALL printMatrix(P, P, V, LDV )
*      print *, "V = Matrix; "
*      print *, "Matrix = []; % Q = "
*      CALL printMatrix(N, N, Q, LDQ )
*      print *, "Q = Matrix; "
*      print *, "Matrix = []; % A_afterSVP = "
*      CALL printMatrix(M, N, A, LDA )
*      print *, "A_afterSVP = Matrix; "
*      print *, "Matrix = []; % B_afterSVP = "
*      CALL printMatrix(P, N, B, LDB )
*      print *, "B_afterSVP = Matrix; "
*      pause
*
*
*     2. Compute the QR factorization
*
      S_QRJ = SECOND()
*
      CALL SGGQRJ( T, L, A( K+1, N-L+1 ), LDA, B( 1, N-L+1 ), LDB,
     $             WORK( IQ1 ), LDQ1, WORK( IQ2 ), LDQ2, WORK( IWK ),
     $             LWK, INFO )
*
      S_QRJ = SECOND() - S_QRJ
*
*      ---- CHECK Result ----
*      print *, "Matrix = []; % A_afterQRJ = "
*      CALL printMatrix(T, L, A( K+1, N-L+1 ), LDA )
*      print *, "A_afterQRJ = Matrix; "
*      print *, "Matrix = []; % B_afterQRJ = "
*      CALL printMatrix(L, L, B( 1, N-L+1 ), LDB )
*      print *, "B_afterQRJ = Matrix; "
*
*      print *, "Matrix = []; % WORKQ1 = "
*      CALL printMatrix(T, L, WORK( IQ1 ), LDQ1 )
*      print *, "WORKQ1 = Matrix; "
*      print *, "Matrix = []; % WORKQ2 = "
*      CALL printMatrix(L, L, WORK( IQ2 ), LDQ2 )
*      print *, "WORKQ2 = Matrix; "
*      pause
*
*     3. Compute the CSD of WORKU and WORKV
*
      CALL SLASET( 'All', N,  1, ZERO, ZERO, ALPHA, 1 )
      CALL SLASET( 'All', N,  1, ZERO, ZERO, BETA,  1 )
      IF( K .GT. 0 ) THEN
         CALL SLASET( 'All', K,  1, ONE,  ONE,  ALPHA, 1 )
*         CALL SLASET( 'All', K,  1, ZERO, ZERO, BETA,  1 )
      END IF
*      CALL SLASET( 'All', G-K,   1, ZERO, ZERO, ALPHA( K+1 ), 1 )
      CALL SLASET( 'All', G-K,   1, ONE, ONE,   BETA( K+1 ),  1 )
*
*
*
      S_CSD = SECOND()
*
      CALL SORCSD( WANTALL, T, L, L, WORK( IQ1 ), LDQ1,
     $             WORK( IQ2 ), LDQ2, ALPHA( K+1 ), BETA( K+1 ),
     $             WORK( IU ),  LDWU,  WORK( IV ), LDWV, WORK( IZ ),
     $             LDWZ, WORK( IWK ), LWK, INFO)
*
      S_CSD = SECOND() - S_CSD

*
*     4. Update U
*
      S_POS = SECOND()
*
      CALL SGEMM( 'N', 'N', M, T, T, ONE, U( 1, K+1 ), LDU,
     $            WORK( IU ), LDWU, ZERO, WORK( IQ1 ), LDQ1 )
      CALL SLACPY( 'All', M, T, WORK( IQ1 ), LDQ1, U( 1, K+1 ), LDU )
*
*      ---- CHECK Result ----
*      print *, "Matrix = []; % U = "
*      CALL printMatrix(M, M, U, LDU )
*      print *, "UpdatedU = Matrix; "
*      pause
*
*     5. Update V
*
      CALL SGEMM( 'N', 'N', P, L, L, ONE, V( 1, 1 ), LDV,
     $              WORK( IV ), LDWV, ZERO, WORK( IQ2 ), LDQ2 )
      CALL SLACPY( 'All', P, L, WORK( IQ2 ), LDQ2, V( 1, 1 ), LDV )
*
*      IF( K+L .LE. P ) THEN
*         DO I = 1, K
*            CALL SCOPY( P, V( 1, L+I ), 1, V( 1, I ), 1 )
*         END DO
*         CALL SLACPY( 'All', P, L, WORK( IQ2 ), LDQ2,
*     $                V( 1, K+1 ), LDV )
*      ELSE
*         print *, "here!!!!! p = ", P,",K=", K, ", L = ", L
*
*         DO I = 1, P-L
*            CALL SCOPY( P, V( 1, L+I ), 1, V( 1, I ), 1 )
*         END DO
*         CALL SLACPY( 'All', P, L, WORK( IQ2 ), LDQ2,
*     $                V( 1, P-L+1 ), LDV )
*      END IF
*
*      ---- CHECK Result ----
*      print *, "V = ..."
*      CALL ptMatrix(P, P, V, LDV )
*
*     6. Multiply WORKZ' R1, store result in WORKQ2
*

      IF( M-K-L .GE. 0 ) THEN
         CALL SGEMM( 'N', 'N', L, L, L, ONE, WORK( IZ ), LDWZ,
     $                A( K+1, N-L+1 ), LDA, ZERO, WORK( IQ2 ), LDQ2 )
      ELSE
*         print *, "%Multiply WORKZ"
         CALL SLASET( 'All', L,  L, ZERO, ZERO, WORK( IQ2 ), LDQ2 )
         CALL SGEMM( 'N', 'N', L, L, M-K, ONE, WORK( IZ ), LDWZ,
     $               A( K+1, N-L+1 ), LDA, ZERO, WORK( IQ2 ), LDQ2 )
         CALL SGEMM( 'N', 'N', L, L-M+K, L-M+K, ONE,
     $               WORK( IZ+(M-K)*LDWZ ),LDWZ, B( 1, N+M-K-L+1 ),
     $               LDB, ONE, WORK( IQ2+(M-K)*LDQ2 ), LDQ2 )


*         CALL SGEMM( 'N', 'N', L, L, M-K, ONE, WORKZ, LDWZ,
*     $               A( K+1, N-L+1 ), LDA, ONE, WORKQ2, LDQ2 )
*
*         CALL SGEMM( 'N', 'N', L, L-M+K, L-M+K, ONE, WORKZ(1, M-K+1),
*     $               LDWZ, B( 1, N-L+M-K+1 ), LDB, ZERO,
*     $               WORKQ2(1, M+K+1), LDQ2 )
      END IF
*
*      ---- CHECK Result ----
*      print *, "Matrix = []; % WORKZ2 = "
*      CALL printMatrix(L, L, WORK( IZ ), LDWZ )
*      print *, "WORKZ2 = Matrix; "
*      print *, "Matrix = []; % A = "
*      CALL printMatrix(M, N, A, LDA )
*      print *, "A_afterMult2 = Matrix; "
*      print *, "Matrix = []; % B = "
*      CALL printMatrix(P, N, B, LDB )
*      print *, "B_afterMult2 = Matrix; "
*      print *, "Matrix = []; "
*      CALL printMatrix(L, L, WORK( IQ2 ), LDQ2 )
*      print *, "WORKQ2_afterMul2 = Matrix; "
*      pause
*
*     7. RQ factorization on result in WORKQ2
*
      CALL SGERQF( L, L, WORK( IQ2 ), LDQ2, WORK( IWK ),
     $             WORK( IWK+(L+1) ), LWK-L, INFO )
*      print *, "Matrix = []; %WORKQ2 "
*      CALL printMatrix(L, L, WORKQ2, LDWQ2 )
*      print *, "WORKQ2_afterRQ = Matrix; "
*      pause

*
*     8. Update R
*
      IF( M-K-L .GE. 0 ) THEN
         CALL SLACPY( 'Upper', L, L, WORK( IQ2 ), LDQ2,
     $                A( K+1, N-L+1 ),  LDA )
      ELSE
*         print *, "%Update R"
         CALL SLACPY( 'Upper', M-K, L, WORK( IQ2 ), LDQ2,
     $                A( K+1, N-L+1 ), LDA )
         I = IQ2+LDQ2*(M-K)+(M-K)
         CALL SLACPY( 'Upper', K+L-M, K+L-M, WORK( I ), LDQ2,
     $               B( 1, N+M-K-L+1 ), LDB )

*         CALL SLACPY( 'Upper', K+L-M, K+L-M, WORKQ2(M+K+1, M+K+1),
*     $               LDWQ2, B( 1,  N-L+M-K+1 ), LDB )
*
      END IF
      CALL SORMR2( 'R', 'T', K, L, L, WORK( IQ2 ), LDQ2, WORK( IWK ),
     $             A( 1, N-L+1 ), LDA, WORK( IWK+L+1 ), INFO )

*
*      ---- CHECK Result ----
*      print *, "Matrix = []; % WORKQ2 = "
*      CALL printMatrix(L, L, WORK( IQ2 ), LDQ2 )
*      print *, "WORKQ2_2 = Matrix; "
*      print *, "Matrix = []; % Rb = "
*      CALL printMatrix(P, N, B, LDB )
*      print *, "Rbf_2 = Matrix; "
*      pause
*
*     9. Update Q
*
      CALL SORMR2( 'R', 'T', N, L, L, WORK( IQ2 ), LDQ2, WORK( IWK ),
     $             Q( 1, N-L+1 ), LDQ, WORK( IWK+L+1 ), INFO )
*     $             Q( N-L+1, 1 ), LDQ, WORK( L+1 ), INFO )
*
*
      S_POS = SECOND() - S_POS
      S_ALL = SECOND() - S_ALL
*
      PRINT *, "in SGGQSV.f, TIME: ", K, ",", L
      WRITE( *, FMT = 9990 ) M, P, N,
     $     S_PRE, 100*S_PRE / S_ALL,
     $     S_QRJ, 100*S_QRJ / S_ALL,
     $     S_CSD, 100*S_CSD / S_ALL,
     $     S_POS, 100*S_POS / S_ALL , S_ALL

 9990  FORMAT(  I4, ' & ', I4, ' & ', I4,
     $          ' & ',  F8.5, ' & ',  F8.2, "\%",
     $          ' & ',  F8.5, ' & ',  F5.2, "\%",
     $          ' & ',  F8.5, ' & ',  F5.2, "\%",
     $          ' & ',  F8.5, ' & ',  F5.2, ' & ', F8.5, '\\')
*
      RETURN
*
*     End of SGGQSV
*
      END
