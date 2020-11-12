      SUBROUTINE SORCS1 ( JOB, M, P, L, Q1, LDQ1, Q2,
     $                    LDQ2, ALPHA, BETA, U, LDU, V, LDV, ZT,
     $                    LDZ, WORK, LWORK, INFO )
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
* SORCS1 computes the cosine-sine decomposition (CSD) of an
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
* D1'D1 + D2'D2 = I.   For details, refer to SORCSD.f.
*
* This is an internal subroutine called by SORCSD to handle the 
* case where assumes M > P.  For other case, see SORCS2.f.
*
* Arguments
* =========
*
* JOB    (input) CHARACTER1
*        = 'Y':  All three orthogonal matrices, U, V, Z are computed;
*        = 'N':  None of orthogonal matrices, U, V, Z are not computed.
*
* M       (input) INTEGER
*         The number of rows of the matrix Q1.  M >= 0.
*
* P       (input) INTEGER
*         The number of rows of the matrix Q2.  P >= 0.
*
* L       (input) INTEGER
*         The number of columns of the matrices Q1 and Q2.
*         L >= 0, and M+P >= L.
*
* Q1      (input/output) REAL array, dimension (LDQ1,L)
*         On entry, the M-by-L matrix Q1.
*         On exit, Q1 is destroyed.
*
* LDQ1    (input) INTEGER
*         The leading dimension of the array Q1. LDQ1 >= max(1,M).
*
* Q2      (input/output) REAL array, dimension (LDQ2,L)
*         On entry, the P-by-L matrix Q2.
*         On exit, Q2 is destroyed.
*
* LDQ2    (input) INTEGER
*         The leading dimension of the array Q2. LDQ1 >= max(1,P).
*
* ALPHA   (output) REAL array, dimension (L)
* BETA    (output) REAL array, dimension (L)
*     On exit, ALPHA and BETA contain the cosine-sine value pairs
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
*     (3) if M < L and P < L
*           ALPHA(1:L-P) = 1, ALPHA(L-P+1:M) = C, ALPHA(M+1, L) = 0
*           BETA (1:L-P) = 0, BETA (L-P+1:M) = S, BETA (M+1, L) = 1
*      
*       ALPHA are stored in non-increasing order.
*       BETA  are stored in non-decreasing order.
*
* U       (output) REAL array, dimension (LDU,M)
*         If JOB = 'Y', U contains the M-by-M orthogonal matrix U.
*         If JOB = 'N', U is not referenced.
*
* LDU     (input) INTEGER
*         The leading dimension of the array U. LDU >= max(1,M) if
*         JOB = 'Y'; LDU >= 1 otherwise.
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
*         If JOB = 'Y', ZT contains the L-by-L orthogonal matrix Z'.
*         If JOB = 'N', ZT is not referenced.
*
* LDZ     (input) INTEGER
*         The leading dimension of the array ZT.
*         LDZ >= max(1,L) if JOB = 'Y'; LDZ >= 1 otherwise.
*
* WORK    (workspace) REAL array, dimension (LWORK)
*
* LWORK   (input) INTEGEER
*         The dimension of the array WORK.
*         LWORK >= MAX( 3 * MIN( M, L ) + MAX( M, L ), 5 * MIN( M, L ) )
*                 + 2 * MIN( P, L ) + L**2 
*
* INFO    (output) INTEGER
*         = 0:  successful exit
*         < 0:  if INFO = -i, the i-th argument had an illegal value.
*
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
*
*     .. Parameters ..
      REAL               ZERO, ONE 
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            WANT
      INTEGER            I, J, R, q_1, q_2, LW, IW, k_p, p_k 

      REAL               TOL, S_SVD1, S_SVD2, S_QR1, S_QR2, S_ALL
*     ..
*     .. External Functions
      LOGICAL            LSAME
      REAL               SECOND 
      EXTERNAL           LSAME, SECOND
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGESVD, SGEQRF, SORM2R
*     ..
*     .. Intrinsic Functions
      INTRINSIC          MAX, MIN, SQRT
*     ..     
*     .. Executable Statements ..
*     Test the input parameters
*
      WANT = LSAME( JOB, 'Y')
*
      INFO = 0
      IF( .NOT. (WANT .OR. LSAME( JOB, 'N') ) ) THEN
        INFO = -1 
      ELSE IF( M.LT.0 ) THEN
        INFO = -2
      ELSE IF( P.LT.0 .OR. M.LE.P ) THEN
        INFO = -3 
      ELSE IF( L.LT.0 ) THEN
        INFO = -4
      ELSE IF( L.LT.0 .OR. ( ( M + P ) .LT. L ) ) THEN
        INFO = -4 
      ELSE IF( LDQ1.LT.MAX( 1, M ) ) THEN
        INFO = -6
      ELSE IF( LDQ2.LT.MAX( 1, P ) ) THEN
        INFO = -8
      ELSE IF( LDU.LT.1 .OR. ( WANT .AND. LDU.LT.M ) ) THEN
        INFO = -12
      ELSE IF( LDV.LT.1 .OR. ( WANT .AND. LDV.LT.P ) ) THEN
        INFO = -14
      ELSE IF( LDZ.LT.1 .OR. ( WANT .AND. LDZ.LT.L ) ) THEN
        INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'SORCS1', -INFO )
        RETURN
      END IF
*
      S_SVD1 = ZERO
      S_SVD2 = ZERO
      S_QR1 = ZERO
      S_QR2 = ZERO
      S_ALL = ZERO
*
*     Test the input parameters
*
      WANT = LSAME( JOB, 'Y')
*
*      S5 = SECOND()
*
*     Set up the threshold for separating the "large" and "small" 
*     singular values.
*      
      TOL = SQRT(0.5)
*
*     Initialize variables
*
      q_1 = MIN( M, L )
      q_2 = MIN( P, L ) 
      p_k = MAX(P-L, 0)
      k_p = MAX(L-P, 0)
      S_ALL = SECOND() 
*   
*--------
*     STEP 1. SVD OF Q1 
*             Q1 = U ALPHA Z', ALPHA is ordered non-increasingly.
*             Free space after this step: Q1, BETA, V, WORK.
*--------
*
      IF( WANT ) THEN
*  
         S_SVD1 = SECOND()
*
         CALL SGESVD( 'A', 'A', M, L, Q1, LDQ1, ALPHA, U, LDU, ZT, 
     $             LDZ, WORK, LWORK, INFO )   
*
         S_SVD1 = SECOND() - S_SVD1
*          
      ELSE
*
*         S_SVD1 = SECOND()
*
         CALL SGESVD( 'N', 'A', M, L, Q1, LDQ1, ALPHA, U, LDU,  
     $             WORK, L, WORK( L*L + 1 ), LWORK-L*L, INFO )   
*
*         S_SVD1 = SECOND() - S_SVD1
*          
      END IF
*
*     1.2 Place ALPHA in BETA while swapping its ordering
*         Fill in zeros and ones in proper places. 
*         May use SLACPY later      
*
      DO 10 I = 1, L - q_2
        ALPHA( I ) = ONE
        BETA ( I ) = ZERO
   10 CONTINUE
*
      DO 20 I = q_1 + 1, L
        ALPHA( I ) = ZERO
        BETA ( I ) = ONE
   20 CONTINUE
*	  
*============= PRINT RESULT AFTER STEP 1 ==================  
*
*      % after STEP 1:  SVD of Q1"        
*      print *, "ZT = ..."
*      CALL PTMATRIX( L, L, ZT, LDZ )
*      print *, "U = ..."
*      CALL PTMATRIX( M, M, U, LDU )
*      print *, " ALPHA = ..."
*      CALL PTMATRIX( L, 1, ALPHA, L )
*      check: U'*Q1*ZT' - diag(beta)      
*===========================================================
*
*--------     
*    STEP 2. 
*     Find R such that 
*      ALPHA(1:R) are the "large" singular values
*      BETA(R+1:L) are the "small" singular values.
*     Free space after this step: Q2, V, WORL(q_2+q_2+1:end)
*--------
*
      R = 0
      IF( ALPHA( q_1 ) .GT. TOL ) THEN
        R = q_1 - k_p
      ELSE
        J = k_p + 1
        DO 30 I = J, q_1
          IF( ALPHA( I ) .LT. TOL ) THEN
            R = I-J
            GOTO 40
          END IF
  30    CONTINUE
      END IF
  40  CONTINUE
      IW = R+k_p
      LW = ( R+k_p )**2
*
*  
*=============== PRINT RESULT AFTER STEP 2 =================
*      print *, "% after STEP 4..., R = ", R, " IW = ", IW
*===========================================================      
*	  
*--------
*     STEP 3. Matrix-matrix multiply, T = Q2 (ZT)'
*             Result stored in Q1
*             Free space after this step: Q2, V, WORK.
*
      IF( WANT ) THEN
        CALL SGEMM( 'N', 'T', P, L, L, ONE, Q2, LDQ2, ZT, LDZ, ZERO, 
     $             Q1, LDQ1 )
      ELSE
        CALL SGEMM( 'N', 'T', P, L, L, ONE, Q2, LDQ2, WORK, L, ZERO, 
     $             Q1, LDQ1 )
      END IF
*
*=============== PRINT RESULT AFTER STEP 3 =================
*
*      print *, " % after STEP 3, T = Q2 * Z ..."
*      print *, "T = ..."
*      CALL PTMATRIX( P, L, Q1, LDQ1 )
*      check:  norm(Q1*Zf-Tf)              
*
*===========================================================
*
*--------
*     STEP 4. QR on T: T = H L
*             Sort L so that its diagonal entries are non-negative.
*             Result H and L stored in Q1,
*             TAU in WORK(1:q_2), sign change info stored in WORK(q2+1:q2+q_2)
*             Free space after this step: Q2, V, WORK(q_2+q_2+1:end)
*--------
*
*     4.1 QR Decomposition of T:  T = H L, 
*         store result in Q2, use WORK(q_2+1:end) as work space, TAU in WORL(1:q_2)
*
      S_QR1 = SECOND()
*
*
      CALL SGEQLF( P, L, Q1, LDQ1, WORK, WORK( q_2+1 ), 
     $             LWORK - q_2, INFO )
*      print *, " Q1 after QL = "
*      CALL PTMATRIX( 1, q_2, WORK, 1 )       
      S_QR1 = SECOND() - S_QR1
*
*     4.2 Make all diagonal of L positive, store info in WORK(q_2+1:q_2+q_2)
      DO 50 I = 1, q_2
        IF( Q1( I+p_k, I+k_p ) .LT. ZERO ) THEN
          CALL SSCAL(I+k_p, -ONE, Q1( I+p_k, 1 ), LDQ1)
          WORK( I+q_2 ) = -ONE
        ELSE
          WORK( I+q_2 ) = ONE 
        END IF
   50   CONTINUE
*  
*=============== PRINT RESULT AFTER STEP 4=================
*
*       print *, " % after STEP 4: QL of T: Q1 = Q1 ZT, Q1(i,i) > 0"
*       print *, " % p_k = ", p_k, " k_p = ", k_p , "P = ", P, "L = ", L
*       print *, " Q1 = ..."
*       CALL PTMATRIX( P, L, Q1, LDQ1 )        
*       print *, " FsignInfo = ... "
*       CALL PTMATRIX( q_2, 1, WORK(q_2+1), 1 )
*      
*===========================================================
*      
*--------
*     STEP 5. Take the first R x (R+k_p) corner of L out, L33, place in Q2
*             Then do the SVD of L33:  V1 S1 Z1' = R33
*             Store V1 in U(q_1+1-T:q_1, q_1+1-T:q_1), 
*             Store Z1' in WORK(q_2+q_2+1, q_2+q_2+(r+k_p)**2).
*             Free space after this step: WORK(q_2+q_2+(r+k_p)**2+1:end)
*--------
*
*     5.0 R < 1, don't need to compute further
*     Formulate BETA and V then return
*      
      IF (R .LT. 1 ) THEN
        DO 60 I = 1 + k_p, q_1
           BETA( I ) = Q1( I+p_k-k_p, I )
  60    CONTINUE
*        
        IF ( WANT ) THEN
*           
          CALL SLASET( 'Full', P, P, ZERO, ZERO, V, LDV )
          CALL SLACPY( 'UPPER', P, q_2, Q1( 1, k_p+1 ), LDQ1, V, LDV )
          CALL SORGQL( P, P, q_2, V( 1, 1 ), LDV, WORK, 
     $                 WORK( q_2+q_2+1 ), LWORK-q_2-q_2, INFO )
*
          DO 70 I = 1, q_2
            IF( WORK( q_2 + I ) .LT. ZERO ) THEN
              CALL SSCAL( P, -ONE, V( 1, I ), 1 )
            END IF
  70      CONTINUE
            print *, "early return "
*           To return
           GOTO 900
        END IF
      END IF  
*
*     5.1 Take out the last R x IW lower triangular of L, L33, put it in Q2
*
*     
      J = p_k - k_p
      DO I = IW+1, q_1 
         BETA( I ) = Q1( I + J, I )
      END DO
 
      CALL SLASET( 'FULL', R, IW, ZERO, ZERO, Q2, LDQ2 )
*      
      CALL SLACPY('FULL ', R, k_p, Q1, LDQ1, Q2, LDQ2 )
      CALL SLACPY('LOWER', R, R,   Q1( p_k+1, k_p+1 ), LDQ1, 
     $            Q2( 1, k_p+1 ), LDQ2) 
*       print *, " L11 = ..."
*       CALL PTMATRIX( R, IW, Q2, LDQ2 )      
     
      IF( WANT ) THEN
*
*       5.2 Initialize partial of V to the identity
*
        IF( p_k .GT. 0 ) THEN 
           CALL SLASET( 'Full', p_k, L,   ZERO, ZERO, V, LDV )
           CALL SLASET( 'Full', p_k, p_k, ZERO, ONE,  V( 1, L+1 ), LDV )
           CALL SLASET( 'Full', R  , P-R, ZERO, ZERO, V( p_k+1, R+1 ), 
     $                  LDV )
           CALL SLASET( 'Full', L-R, R,   ZERO, ZERO, V( p_k+R+1, 1 ), 
     $                  LDV )
           CALL SLASET( 'Full', L-R, L-R, ZERO, ONE,  V( p_k+R+1, R+1 ),
     $                  LDV )
           CALL SLASET( 'Full', L-R, p_k, ZERO, ZERO, V( p_k+R+1, L+1 ),
     $                  LDV )
        ELSE
           CALL SLASET( 'Full', R, P-R,   ZERO, ZERO, V( 1, R+1 ), LDV )
           CALL SLASET( 'Full', P-R, R,   ZERO, ZERO, V( R+1, 1 ), LDV )
           CALL SLASET( 'Full', P-R, P-R, ZERO, ONE, V( R+1, R+1), LDV )
        END IF
*
*       5.3 Compute the SVD of L33:  V1 S1 Z1' = L33
*         store V1 in portion of V, store Z1' in WORK start at  WORK( q_2+q2+1 )
*
        S_SVD2 = SECOND()
        CALL SGESVD( 'A', 'A', R, IW, Q2, LDQ2, BETA( k_p + 1 ), 
     $       V( p_k+1, 1 ), LDV, WORK(q_2+q_2+1), IW, 
     $       WORK(q_2+q_2+LW+1), LWORK-LW-q_2-q_2, INFO )
*
        S_SVD2 = SECOND() - S_SVD2
*
*     5.3 Swap the columns of V1 and rows of Z1
*      
*       Swap columns of V1
*
      DO 80 I = 1, R/2
          CALL SSWAP( R, V( p_k+1, I ), 1, V( p_k+1, R-I+1 ), 1 )
   80   CONTINUE
*
*       SWAP rows of ZT
*
        J = q_2 + q_2
        DO 90 I = 1, IW/2
          CALL SSWAP( IW, WORK( J+I ), IW, WORK( J+IW-I+1 ), IW )
   90   CONTINUE      
*          
      ELSE
*
        CALL SGESVD( 'N', 'N', R, R+k_p, Q2, LDQ2, BETA( k_p + 1 ), 
     $       V, LDV, U, LDU, WORK, LWORK, INFO )
*
      END IF
*
*     SWAP the elements in BETA(k_p+1:k_p+R)
      DO I = 1, R/2
        CALL SSWAP( 1, BETA( k_p + I ), 1, BETA( R + k_p + 1 - I ), 1 ) 
      END DO
*
*
*=============== PRINT RESULT AFTER STEP 5 =================
*
*       print *, " % after STEP 5: SVD of L11: V1*S1*Z1=L11 and swap "
*       print *, " V1 = ..."
*       CALL PTMATRIX( R, R, V(p_k+1, 1), LDV )
*       print *, " Z1 = ..."
*       CALL PTMATRIX( IW, IW, WORK(q_2+q_2+1), IW ) 
*       print *, " BETA = ..."
*       CALL PTMATRIX( 1, R, BETA( k_p + 1 ), 1 )
*       print *, " LargeV = ..."
*       CALL PTMATRIX( P, P, V, LDV )
*
*       print *, "ALPHA = ..."
*       CALL PTMATRIX( 1, L, ALPHA, 1 )
*       print *, "BETA = ..."       
*       CALL PTMATRIX( 1, L, BETA , 1 )       
*        
*===========================================================
*
*--------
*     STEP 6. Optional: formulate V
*     Free space after this step: Q1
*--------
*      
      IF( WANT ) THEN
*
*	  6.1. Change the signs in the last q_2 rows of V
*              accordingly to WORK(1:q_2)
*
        DO 100 I = 1, q_2
*          print *, "WORK = ", WORK( I + q_2 )
          IF( WORK( I + q_2 ) .LT. ZERO ) THEN
             CALL SSCAL(P, -ONE, V(I+p_k, 1), LDV)
          END IF
  100   CONTINUE
*		     
*     6.2: V = V V1 where H is the Householder reflector from 3.1.
*
*        CALL SORGQL(P, q_2, q_2, Q1( 1, k_p+1 ), LDQ1, WORK, 
*     $              WORK(LW+q_2+q_2+1), P, INFO )
*        CALL PTMATRIX( P, q_2, Q1( 1, k_p+1 ), LDQ1 )              
        CALL SORMQL( 'L','N', P, P, q_2, Q1( 1, k_p+1 ), LDQ1, WORK, 
     $              V, LDV, WORK(LW+q_2+q_2+1), P, INFO)
*	 
*  
*=============== PRINT RESULT AFTER STEP 6 =================
*
*        print *, "% after STEP 6: formulate V ..."
*        print *, " V = ... "
*        CALL PTMATRIX( P, P, V, LDV )
*===========================================================
*        
*
*--------
*     STEP 7. Optional: formulate ZT
*     Free space after this step: Q2, WORK(1:end)
*--------
*
*     7.1. Q2 = Z1'*Z(1:L, 1:R+k_p)', use Q1 as work space
*          note:  Q2 has sufficient space since M >= R-k_p.        
*
*
        I = MIN( IW, M )
        J = 1
  110   CONTINUE        
        
        CALL SLASET( 'Full', I, L, ZERO, ZERO, Q1, LDQ1 )
        CALL SGEMM( 'N', 'N', I, L, I, ONE, WORK(q_2+q_2+1), I,
     $          ZT( J, 1 ), LDZ, ZERO, Q1, LDQ1 )
  
*     7.2. Move the resulted part from Q1 to ZT, 
*
      CALL SLACPY('Full', I, L, Q1(1, 1), LDQ1, ZT( J, 1 ), LDZ )
*  
*=============== PRINT RESULT AFTER STEP 8 =================
*
*        print *, "% after STEP 7: formulate ZT ... "
*        print *, " ZT = ..." 
*        CALL PTMATRIX(L, L, Z, LDZ)        
*===========================================================
*
*--------
*     STEP 8. Optional: Formulate U
*     Free space after this step: Q1,Q2, WORK(1:end).
*--------
*
*     Want D*Q2, but have Q2'.  So D*Q2 = (D')*(Q2')'=(Q2' * D)' 
*     X * D = X(:, i) * D(i,i) for each i in columns of X and D
*
*     8.1. get the first IW values of ALPHA and multiply to Z1'
*        call it W, stored in WORK
* 
        DO 120  J = 1, IW
         CALL SSCAL( IW, ALPHA( J ), WORK( q_2+q_2+1+(J-1)*IW ), 1)
  120   CONTINUE
*
*     8.2. Take the transpose of W
*          store result in Q1
*
      CALL SLASET( 'FULL', IW, IW, ZERO, ZERO, Q1, LDQ1 )
      DO 140 I = 1,IW
        DO 130 J = 1,IW
          Q1(I, J) = WORK( q_2+q_2+(I-1)*IW+J )
*                  = WORK( J, I )          
  130   CONTINUE
  140 CONTINUE
*      print *, "W = .."
*	  call ptmatrix(IW, IW, Q1(1,1), LDQ1)
*      
*     8.3. Do QR on W:  W = Uw Rw, stored in Q1, TAU in WORK( 1:IW )
*
      S_QR2 = SECOND()
      CALL SGEQRF( IW, IW, Q1, LDQ1, WORK, WORK( IW + 1 ), 
     $             LWORK-IW, INFO )
*
      S_QR2 = SECOND() - S_QR2
*
*     8.4. Make sure the diagonal entries are positive, 
*          Save additional sign info in WORK( IW+1, IW+IW )
*      
        DO 150 I = 1, IW
        IF( Q1( I, I ) .LT. ZERO ) THEN
           CALL SSCAL( IW-I+1, -ONE, Q1( I, I ), LDQ1 )
           WORK( IW+I ) = -ONE
        ELSE
          WORK( IW+I ) = ONE
        END IF
  150   CONTINUE
*
*     8.5. Compute U = U * Qw, 
*          where Qw is the Houlseholder reflectors.
*
        CALL SORM2R( 'R', 'N', M, IW, IW, Q1, LDQ1, WORK, U,
     $               LDU, WORK( IW*2+1 ), INFO )        
*
*     8.6. Change the signs of U(1:m 1:IW ) according to 
*          WORK( IW+1:IW+IW )
*
      DO 170 J = 1, IW
        IF( WORK( IW + J ) .LT. ZERO ) THEN
           CALL SSCAL( M, -ONE, U( 1, J ), 1 )
*         DO 160 I = 1, M
*           U( I, J ) = -U( I, J )
* 160     CONTINUE
        END IF
  170 CONTINUE

      END IF
*  
*=============== PRINT RESULT AFTER STEP 8 =================
*
*      print *, "% after STEP 8: formulate U ..."
*      print *, "U = ..."
*      CALL PTMATRIX( M, M, U, LDU )
*===========================================================
*
*     
  900 CONTINUE
       S_ALL = SECOND() - S_ALL
      PRINT *, "in SORCS1.f, tIME: "
      WRITE( *, FMT = 9990 ) M, P, L, S_SVD1, S_SVD2, 
     $     (100*( S_SVD1 + S_SVD2 ) / S_ALL),  S_QR1, S_QR2, 
     $     (100*(S_QR1 + S_QR2)/S_ALL), S_ALL

 9990  FORMAT(  I4, ' & ', I4, ' & ', I4, 
     $          ' & ',  F8.5, ' & ',  F8.5, 
     $          ' & ',  F8.2, ' & ',  F8.5, 
     $          ' & ',  F8.5, ' & ',  F8.2, ' & ', F8.5, ' \\') 
*
      RETURN
*
*     End of SORCS1
*
      END      
