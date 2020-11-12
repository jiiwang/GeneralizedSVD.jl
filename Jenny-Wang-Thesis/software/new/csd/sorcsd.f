      SUBROUTINE SORCSD ( JOB, M, P, L, Q1, LDQ1, Q2,
     $                    LDQ2, ALPHA, BETA, U, LDU, V, LDV, ZT,
     $                    LDZ,  WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK driver routine (version 0.1) --
*     Univ. of California, Davis
*     September, 2005
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            INFO, LDQ1, LDQ2, LDU, LDV, LDZ, M, P,
     $                   L, LWORK
*     ..
*     .. Array Arguments ..
      REAL               Q1( LDQ1, * ), ALPHA( * ), Q2( LDQ2, * ),
     $                   BETA( * ), U( LDU, * ), V( LDV, * ),
     $                   ZT( LDZ, * ), WORK( LWORK )
*     ..
*
* PURPOSE:
* ========
*
* SORCSD computes the cosine-sine decomposition (CSD) of an
* (M+P)-by-L orthogonal matrix Q = (Q1' Q2')':
*
*       U' * Q1 = D1 Z',   V' * Q2 = D2 Z'
*
*
* where U, V and Z are orthogonal matrices of size M-by-M, P-by-P,
* and L-by-L, respectively.  X' denote the transpose of X.
* ( NOTE: the actual return value is Z' instead of Z )
*
* D1 and D2 are M-by-L and P-by-L "diagonal" matrices and
* D1'D1 + D2'D2 = I.  D1 and D2 have one of the following four
* structures:
*
* (1) If M >= L and P >= L,
*
* D1 =        L
*         L ( C )
*       M-L ( 0 )
*
* D2 =        L
*         L ( S )
*       P-L ( 0 )
*
* where C = diag(ALPHA(1), ... , ALPHA(L) ),
*       S = diag(BETA(1), ..., BETA(L) ),
*       C**2 + S**2 = I.
*
* (2) If M >= L and P < L,
*
* D1 =          L-P  P
*        L-P  (  I   0 )
*          P  (  0   C )
*        M-L  (  0   0 )
*
* D2 =         L-P   P
*          P  ( 0    S )
*
* where  C = diag(ALPHA(L-P+1), ... , ALPHA(L) ),
*        S = diag(BETA(L-P+1), ..., BETA(L) ),
*        C**2 + S**2 = I.
*
* (3) If M <= L and P >= L,
*
* D1 =       M   L-M
*       M  ( C    0  )
*
* D2 =       M    L-M
*       M  ( S     0  )
*     L-M  ( 0     I  )
*     P-L  ( 0     0  )
*
* where C = diag(ALPHA(1), ... , ALPHA(M) ),
*       S = diag(BETA(1), ..., BETA(M) ),
*       C**2 + S**2 = I.
*
* (4) If M <= L and P < L,
*
* D1 =         L-P   M+P-L   L-M
*        L-P (  I      0      0  )
*      M+P-L (  0      C      0  )
*
* D1 =         L-P   M+P-L   L-M
*      M+P-L (  0      S      0  )
*        L-M (  0      0      I  )
*
* where C = diag(ALPHA(L-P+1), ... , ALPHA(M) ),
*       S = diag(BETA(L-P+1), ..., BETA(M) ),
*       C**2 + S**2 = I.
*
*
* Arguments
* =========
*
* JOB     (input) CHARACTER*1
*         = 'Y':  Orthogonal matrix U,V, and Z are computed;
*         = 'N':  Orthogonal matrices are not computed.
*
* M       (input) INTEGER
*         The number of rows of the block of Q, Q1.  M >= 0.
*
* P       (input) INTEGER
*         The number of rows of the block of Q, Q2.  P >= 0.
*
* L       (input) INTEGER
*         The number of columns of Q.
*         L >= 0, and M+P >= L.
*
* Q1      (input/output) REAL array, dimension (LDQ1,*)
*         On entry, the M-by-L matrix Q1.
*         On exit, Q1 is destroyed.
*
* LDQ1    (input) INTEGER
*         The leading dimension of the array Q1. LDQ1 >= max(1,M).
*
* Q2      (input/output) REAL array, dimension (LDQ2,*)
*         On entry, the P-by-L matrix Q2.
*         On exit, Q2 is destroyed.
*
* LDQ2    (input) INTEGER
*         The leading dimension of the array Q2. LDQ1 >= max(1,P).
*
* ALPHA   (output) REAL array, dimension (L)
* BETA    (output) REAL array, dimension (L)
*     On exit, ALPHA and BETA contain the cosine-sine values
*         of Q1 and Q2
*
*     (1) if M >= L and P >= L
*           ALPHA(1:L) = C
*           BETA (1:L) = S
*
*     (2) if M >= L and P < L
*           ALPHA(1:L-P) = 1, ALPHA(L-P+1:L) = C
*           BETA (1:L-P) = 0, BETA (L-P+1:L) = S
*
*     (3) if M <= L and P >= L
*           ALPHA(1:M) = C, ALPHA(M+1:L) = 0
*           BETA (1:M) = S, BETA (M+1:L) = 1
*
*     (4) if M <= L and P < L
*           ALPHA(1:L-P) = 1, ALPHA(L-P+1:M) = C, ALPHA(M+1, L) = 0
*           BETA (1:L-P) = 0, BETA (L-P+1:M) = S, BETA (M+1, L) = 1
*
*       ALPHA are stored in non-increasing order.
*       BETA  are stored in non-decreasing order.
*
* U       (output) REAL array, dimension (LDU,M)
*         If JOB = 'U', U contains the M-by-M orthogonal matrix U.
*         If JOB = 'N', U is not referenced.
*
* LDU     (input) INTEGER
*         The leading dimension of the array U. LDU >= max(1,M)
*         if JOB = 'Y'; LDU >= 1 otherwise.
*
* V       (output) REAL array, dimension (LDV,P)
*         If JOB = 'Y', V contains the P-by-P orthogonal matrix V.
*         If JOB = 'N', V is not referenced.
*
* LDV     (input) INTEGER
*         The leading dimension of the array V. LDV >= max(1,P) if
*         JOB = 'Y'; LDV >= 1 otherwise.
*
* ZT       (output) REAL array, dimension (LDZ,L)
*         If JOB = 'Y', Z contains the L-by-L orthogonal matrix Z'.
*         If JOB = 'N', Z is not referenced.
*
* LDZ     (input) INTEGER
*         The leading dimension of the array ZT.
*         LDZ >= max(1,L) if JOB = 'Y'; LDZ >= 1 otherwise.
*
* WORK    (workspace) REAL array, dimension (LWORK)
*
* LWORK   (input) INTEGEER
*         The dimension of the array WORK.
*         LWORK >= MAX( 3 * MIN( N ,L ) + MAX( N, L ), 5 * MIN( N,L ) )
*                 + 2 * MIN( M, P, L ) + L**2
*
*          where N = MAX( M, P )
*
* INFO    (output) INTEGER
*         = 0:  successful exit
*         < 0:  if INFO = -i, the i-th argument had an illegal value.
*
* Internal Parameters
* ===================
*
* TOL    REAL
*        Thresholds to separate the "large" and "small" singular
*        values.
*        Currently, it is set to SQT(0.5).
*        The size of TOL has the effect of minimizing the backward
*           error of the composition.
*
* ==================================================================
*
*     .. Local Scalars ..
      LOGICAL            WANTALL
*     ..
*     .. External Functions
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SORCS1, SORCS2, XERBLA
*     ..
*     .. Intrinsic Functions
      INTRINSIC          MAX
*     ..
*     Test the input parameters
*
      WANTALL = LSAME( JOB, 'Y')
*
      INFO = 0
      IF( .NOT. (WANTALL .OR. LSAME( JOB, 'N') ) ) THEN
        INFO = -1
      ELSE IF( M.LT.0 ) THEN
        INFO = -2
      ELSE IF( P.LT.0 ) THEN
        INFO = -3
      ELSE IF( L.LT.0 ) THEN
        INFO = -4
      ELSE IF( L.LT.0 .OR. ( ( M + P ) .LT. L ) ) THEN
        INFO = -4
      ELSE IF( LDQ1.LT.MAX( 1, M ) ) THEN
        INFO = -6
      ELSE IF( LDQ2.LT.MAX( 1, P ) ) THEN
        INFO = -8
      ELSE IF( LDU.LT.1 .OR. ( WANTALL .AND. LDU.LT.M ) ) THEN
        INFO = -12
      ELSE IF( LDV.LT.1 .OR. ( WANTALL .AND. LDV.LT.P ) ) THEN
        INFO = -14
      ELSE IF( LDZ.LT.1 .OR. ( WANTALL .AND. LDZ.LT.L ) ) THEN
        INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'SORCSD', -INFO )
        RETURN
      END IF
*
      IF( ( L.EQ.0 ) .OR. ( M.EQ.0 ) .AND. ( P.EQ.0 ) )  RETURN
*
*       PRINT *, "A = ..."
*       CALL PTMATRIX( M, L, Q1, LDQ1 )
*       PRINT *, "B = ..."
*       CALL PTMATRIX( P, L, Q2, LDQ2 )
*
      IF( M .GT. P ) THEN
*
*     -------
*     Case 1: M > P
*     -------
*
       CALL SORCS1 ( JOB, M, P, L, Q1, LDQ1, Q2,
     $               LDQ2, ALPHA, BETA, U, LDU, V, LDV, ZT,
     $               LDZ,  WORK, LWORK, INFO )
*
      ELSE
*
*     -------
*     Case 2: M <= P
*     -------
*
       CALL SORCS2 ( JOB, M, P, L, Q1, LDQ1, Q2,
     $               LDQ2, ALPHA, BETA, U, LDU, V, LDV, ZT,
     $               LDZ,  WORK, LWORK, INFO )
*
      END IF
*      S5 = SECOND() - S5
*      PRINT *, " SORCSD: S1 = , ", S1, ", S2 = ", S2,
*     $         ", S3 = ", S3, ", S4 = ",  S4, " , S5 = ", S5
*      PRINT *, " SORCSD: ( S1 + S2 ) / S5 = ", (S1 + S2)/S5
*      PRINT *, " SORCSD: ( S3 + S4 ) / S5 = ", (S3 + S4)/S5
*
*      print *, " ALPHA = ..."
*      CALL PTMATRIX( L, 1, ALPHA, L )
*      print *, " BETA = ..."
*      CALL PTMATRIX( L, 1, BETA, L )

      RETURN
*
*     End of SORCSD
*
      END
