      SUBROUTINE SORCS2 ( JOB, M, P, L, Q1, LDQ1, Q2,
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
* case where assumes M <= P.  For other case, see SORCS1.f.
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
*     (2) if M < L and P >= L
*           ALPHA(1:M) = C, ALPHA(M+1:L) = 0
*           BETA (1:M) = S, BETA (M+1:L) = 1
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
*         LWORK >= MAX( 3 * MIN( P, L ) + MAX( P, L ), 5 * MIN( P, L ) )
*                 + 2 * MIN( M, L ) + L**2 
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
      INTEGER            I, J, R, q_1, q_2, T, T2, LW, IW

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
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..     
*     .. Executable Statements ..
*     ..     
*     Test the input parameters
*
      WANT = LSAME( JOB, 'Y')
*
      INFO = 0
      IF( .NOT. (WANT .OR. LSAME( JOB, 'N') ) ) THEN
        INFO = -1 
      ELSE IF( M.LT.0 ) THEN
        INFO = -2
      ELSE IF( P.LT.0 .OR. M.GT.P ) THEN
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
        CALL XERBLA( 'SORCS2', -INFO )
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
      S_ALL = SECOND() 
*     
*
*       PRINT *, "======= SORCS2... M = ", M, ", P = ", P, ", L = ", L
*      I = MAX( 3*MIN( P, L ) + MAX( P, L ), 5*MIN( P, L ) ) 
*      print *, "I = ", I, " LWORK = ", LWORK, ", WANT = ", WANT 
*   
*--------
*     STEP 1. SVD OF Q2 and switch ordering of its singular values.
*             Q2 = V BETA Z', BETA is ordered non-decreasingly.
*             Free space after this step: Q2, ALPHA, U, WORK.
*--------
*     1.1 Compute the SVD of Q2: Q2 = V ALPHA Z'
*
      IF( WANT ) THEN
*  
         S_SVD1 = SECOND()
*
         CALL SGESVD( 'A', 'A', P, L, Q2, LDQ2, ALPHA, V, LDV, ZT, 
     $             LDZ, WORK, LWORK, INFO )   
*
         S_SVD1 = SECOND() - S_SVD1
*
*       Swap columns of V
*         
        DO 10 I = 1, q_2/2
          CALL SSWAP( P, V( 1, I ), 1, V( 1, q_2-I+1 ), 1 )
   10   CONTINUE
*
*       SWAP rows of ZT
*
        DO 20 I = 1, L/2
          CALL SSWAP( L, ZT( I, 1 ), LDZ, ZT( L-I+1, 1 ), LDZ )
   20   CONTINUE
*          
      ELSE
*
*         S_SVD1 = SECOND()
*
         CALL SGESVD( 'N', 'A', P, L, Q2, LDQ2, ALPHA, V, LDV,  
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
*      CALL SCOPY( q_1-L+q_2, ALPHA( q_2 ), 1, BETA( L-q_2+1 ), 1 )
      J = 1
      DO 30 I = L - q_2 + 1, q_1 
        BETA( I ) = ALPHA( q_2+1-J ) 
        J = J + 1
   30 CONTINUE   
      DO 40 I = 1, L - q_2
        BETA ( I ) = ZERO
        ALPHA( I ) = ONE
   40 CONTINUE
      DO 50 I = q_1 + 1, L
        BETA ( I ) = ONE
        ALPHA( I ) = ZERO
   50 CONTINUE
*	  
*============= PRINT RESULT AFTER STEP 1 ==================  
*
*      print *, "% after STEP 1: "
*      print *, " ZT = ..."
*      CALL PTMATRIX( L, L, ZT, LDZ )
*      print *, " V = ..." 
*      CALL PTMATRIX( P, P, V, LDV )
*      print *, "BETA = ..."
*      CALL PTMATRIX( L, 1, BETA, L )
*      check: Vf'*Q2*Zf - diag(beta)      
*===========================================================
*
*--------     
*    STEP 2. 
*     Find R such that 
*      BETA(1:R) are the "small" singular values
*      BETA(R+1:q_2) are the "large" singular values.
*     Free space after this step: U, ALPHA, Q1, WORK(q_1+q_1+1:end)
*--------
*
      R = 0
      IF( ( q_2 .GT. q_1 ) .AND. ( BETA (q_1 + 1 ) .LE. TOL ) ) THEN
          R = q_1 
      ELSE
        J = MAX( L -  q_2, 1)
        DO 70 I = J, q_2
          IF( BETA( I ) .LE. TOL )  THEN
            R = I
          ELSE
            GOTO 80
          END IF
  70    CONTINUE
*
  80   CONTINUE
      END IF 
*
      T = q_1 - R
*  
*=============== PRINT RESULT AFTER STEP 2 =================
*      print *, "% after STEP 2:  R = ", R, " T = ", T
*===========================================================      
*	  
*
*--------
*     STEP 3. Matrix-matrix multiply, T = Q1 (Z')'
*             Result stored in Q2
*             Free space after this step: U, ALPHA, Q1, WORK.
*
*      CALL SLASET( 'FULL', L, L, ZERO, ZERO, Q2, LDQ2 )
      IF( WANT ) THEN
        CALL SGEMM( 'N', 'T', M, L, L, ONE, Q1, LDQ1, ZT, LDZ, ZERO, 
     $             Q2, LDQ2 )
      ELSE
        CALL SGEMM( 'N', 'T', M, L, L, ONE, Q1, LDQ1, WORK, L, ZERO, 
     $             Q2, LDQ2 )
      END IF
*
*=============== PRINT RESULT AFTER STEP 3 =================
*
*      print *, "% after STEP 3: "
*      print *, " Q2 = ..."      
*      CALL PTMATRIX( M, L, Q2, LDQ2 )
*      check:  norm(Q1*Zf-Tf)              
*
*===========================================================
*
*--------
*     STEP 4. QR on T: T = H R
*             Sort R so that its diagonal entries are non-negative.
*             Result H and R stored in Q2, 
*             sign change info stored in WORK(1:q_1)
*             Free space after this step: U, ALPHA, Q1, WORK(q_1+q_1:end)
*--------
*
*     4.1 QR Decomposition of T:  T = H R, 
*         store result in Q2, use WORK(q_1+1:end) as work space, TAU in WORK(1:q_1)
*
      S_QR1 = SECOND()
*
*
      CALL SGEQRF( M, L, Q2, LDQ2, WORK, WORK( q_1+1 ), 
     $             LWORK - q_1, INFO )
      
      S_QR1 = SECOND() - S_QR1
*
*     4.2 Make all diagonal of R positive, store info in WORK(L+1:q_1+q_1)
*     May use SSCAL later
*     
      DO 60 I = 1, q_1
        IF( Q2( I, I ) .LT. ZERO ) THEN
*         CALL SSCAL(q_1-I+1, -ONE, Q2(I, I), LDQ2)         
          CALL SSCAL(L-I+1, -ONE, Q2(I, I), LDQ2)                 
          WORK( I+q_1 ) = -ONE
        ELSE
          WORK( I+q_1 ) = ONE 
        END IF
   60 CONTINUE
*  
*=============== PRINT RESULT AFTER STEP 4=================
*       print *, "% after STEP 4: "
*       print *, "Q2 = ..."        
*       CALL PTMATRIX( M, L, Q2, LDQ2 )
*       print *, "FsignInfo = ... "
*       CALL PTMATRIX( q_1, 1, WORK(q_1+1), 1 )
*      
*===========================================================
*      
*--------
*     STEP 5. Take the last T x T corner of R out, R33, place in Q1
*             Then do the SVD of R33:  Ur1 ALPHA Vr1' = R33
*             Store Ur1 in U(q_1+1-T:q_1, q_1+1-T:q_1), 
*             Store Vr1' in leading part of Q1.
*             Free space after this step: WORK(q_1+q_1+1:end)
*--------
*
*     5.0 T < 1, don't need to compute further
*     Formulate ALPHA and U then return
*      
      IF (T .LT. 1 ) THEN
        DO 90 I = 1, L - q_2
           ALPHA( I ) = ONE
   90   CONTINUE
*        I = L - q_2 + 1
*        CALL SCOPY( q_1-I+1, Q2( I, I ), LDQ2, ALPHA( I ), 1 ) 
        DO I = L - q_2 + 1, q_1
           ALPHA( I ) = Q2( I, I )
        END DO
        DO 100 I = q_1 + 1, L
           ALPHA( I ) = ZERO
  100   CONTINUE
*        
        IF ( WANT ) THEN
*           
          CALL SLASET( 'Full', M, M, ZERO, ZERO, U, LDU )
          CALL SLACPY( 'LOWER', M, q_1, Q2, LDQ2, U, LDU )
          CALL SORGQR( M, M, q_1, U, LDU, WORK, WORK( q_1+q_1+1 ), 
     $           LWORK-q_1-q_1, INFO )
*
          DO 110 I = 1, q_1
            IF( WORK( q_1 + I ) .LT. ZERO ) THEN
              CALL SSCAL( M, -ONE, U( 1, I ), 1 )
*              DO J = 1, M
*                U(J, I) = -U(J, I)
*              END DO
            END IF
  110   CONTINUE
*         To return
*          print *, "Early return"
          GOTO 900
        END IF
      END IF  
*
*     5.1 Take out the last TxIW upper triangular of R, R33, put it in Q1      
*
      T2 = q_2 - R
      IW = L - R
      LW = IW*IW
*      
      CALL SLASET( 'FULL', T, IW, ZERO, ZERO, Q1, LDQ1 )
*      
      CALL SLACPY('Upper', T, IW, Q2(R+1, R+1), LDQ2, Q1, LDQ1) 
     
      IF( WANT ) THEN
*
*       5.2 Initialize U to the identity
*
        CALL SLASET( 'Full', M, M, ZERO, ONE, U, LDU )
*
*       5.3 Compute the SVD of R33:  UR1 BETAR1 VR1' R33
*         store Ur1 in U, store Vr1' in Q1 
*
        S_SVD2 = SECOND()
*        CALL SGESVD( 'A', 'A', T, IW, Q1, LDQ1, ALPHA( R + 1 ), 
*     $     U( R + 1, R + 1), LDU, WORKZ, LDWZ, 
*     $     WORK(q_1+q_1+1), LWORK-q_1-q_1, INFO )
        CALL SGESVD( 'A', 'A', T, IW, Q1, LDQ1, ALPHA( R + 1 ), 
     $       U( R + 1, R + 1), LDU, WORK(q_1+q_1+1), IW, 
     $       WORK(q_1+q_1+LW+1), LWORK-LW-q_1-q_1, INFO )
*      
        S_SVD2 = SECOND() - S_SVD2
      ELSE
*
        CALL SGESVD( 'N', 'N', T, IW, Q1, LDQ1, ALPHA( R + 1 ), 
     $       U, LDU, V, LDV, WORK, LWORK, INFO )
*
      END IF
*
      DO 130 I = L - q_2 + 1, R
        ALPHA( I ) = Q2( I, I )
  130 CONTINUE
*  
*=============== PRINT RESULT AFTER STEP 5 =================
*
*      print *, "% after STEP 5:   "
*      print *, "Ur = ... % (SVD result in U( T:q_1, T:q_1) ) "
*      CALL PTMATRIX( T, T, U(R+1,R+1), LDU )
*      print *, "Zr = ..."      
*      CALL PTMATRIX( IW, IW, WORK(q_1+q_1+1), IW )      
*      print *, "Cr = ..."      
*      CALL PTMATRIX( T, 1, ALPHA( R + 1 ), L )
*      print *, " ALPHA = ..."
*      CALL PTMATRIX( L, 1, ALPHA, L )
*      print *, " BETA = ..."
*      CALL PTMATRIX( L, 1, BETA, L )        
*      check : norm(Ur'*R33*Zr''-[diag(Cr(1:m)) zeros(m, k-m)])
*
*============================================================
*
*--------
*     STEP 6. Optional: formulate U
*     Free space after this step: Q2
*--------
*      
      IF( WANT ) THEN
*
*	  6.1. Change the signs of rows of U(1:q_1, 1:M) 
*              accordingly to WORK(1:q_1)
*
        DO 150 I = 1, q_1
*          print *, "WORK = ", WORK( I + q_1 )
          IF( WORK( I + q_1 ) .LT. ZERO ) THEN
             CALL SSCAL(M, -ONE, U(I, 1), LDU)
          END IF
  150   CONTINUE
*		     
*     6.2: U = H Ur1 where H is the Householder reflector from 3.1.
*
        CALL SORM2R( 'L','N', M, M, q_1, Q2, LDQ2, WORK, U, LDU, 
     $              WORK (LW+q_1+q_1+1), INFO)
*	 
*  
*=============== PRINT RESULT AFTER STEP 6 =================
*
*        print *, " % after STEP 6: "
*        print *, " U = ..."        
*        CALL PTMATRIX(M, M, U, LDU )
*===========================================================
*        
*
*--------
*     STEP 7. Optional: formulate ZT
*     Free space after this step: Q2, WORK(1:end)
*--------
*
*     7.1. Q2 = Vr1*ZT(L+1-T:L, :), use Q2 as work space
*
        CALL SLASET( 'Full', IW, L, ZERO, ZERO, Q2, LDQ2 )
*        CALL SGEMM( 'N', 'N', IW, L, IW, ONE, WORKZ, LDWZ,
*     $          Z( R + 1, 1), LDZ, ZERO, 
*     $          Q2, LDQ2 )
        CALL SGEMM( 'N', 'N', IW, L, IW, ONE, WORK(q_1+q_1+1) , IW,
     $          ZT( R + 1, 1), LDZ, ZERO, 
     $          Q2, LDQ2 )
  
*     7.2. Move the resulted part from Q2 to ZT, 
*
      CALL SLACPY('Full', IW, L, Q2(1, 1), LDQ2, ZT( R+1, 1 ), LDZ)   
*  
*=============== PRINT RESULT AFTER STEP 7 =================
*
*        print *, " % after STEP 7:"
*        print *, " ZT = ..."              
*        CALL PTMATRIX(L, L, ZT, LDZ)        
*===========================================================
*
*--------
*     STEP 8. Optional: Formulate V
*     Free space after this step: Q1,Q2, WORK(1:end).
*--------
*
*        IF( 1+R .LE. L-q_2 )  THEN
*           GOTO 900           
*        END IF
*        
*     Want D*Q1, but have Q1'.  So D*Q1 = (D')*(Q1')'=(Q1' * D)' 
*     X * D = X(:, i) * D(i,i) for each i in columns of X and D
*
*     8.1. get the last T2 of beta and multiply to VR1
*        call it W, stored in Q1
* 
        I = q_1+q_1+1
        DO 160  J = 1, T2
*         CALL SSCAL( T2, BETA( R+J ), WORKZ( 1, J ), 1)        
         CALL SSCAL( T2, BETA( R+J ), WORK( q_1+q_1+1+(J-1)*IW ), 1)
*         Start at ZrT(1, L-q_2+1) = WORK( q_1+q_1+(L-q_2+1-1)*IW+1 )
*         CALL SSCAL( T2, BETA( R+J ), WORK( I+(L-q_2+J-1)*IW ), 1)    
  160   CONTINUE
*
*     8.2. Take the transpose of Q1( 1 : T, 1 : T ), 
*          store result in Q2( 1 : T, 1 : T )
*
*     May use SGEMM later
*            
*      CALL SLASET( 'FULL', q_2 - R, L - R, ZERO, ZERO, Q2, LDQ2 )
      CALL SLASET( 'FULL', T2, T2, ZERO, ZERO, Q2, LDQ2 )
      DO 180 I = 1,T2
*        DO J = 1,L - R
        DO 170 J = 1,T2
*          Q2(I, J) = WORKZ(J, I)
          Q2(I, J) = WORK( q_1+q_1+(I-1)*IW+J )
  170   CONTINUE
  180 CONTINUE
*      
*     8.3. Do QR on W:  W = Ur2 T, stored in Q2, TAU in WORK( 1:q_2-R )
*
      S_QR2 = SECOND()
*      CALL SGEQRF( q_2 - R, L - R, Q2, LDQ2, WORK, WORK( q_2 - R + 1 ), 
*     $             LWORK-q_2+R+1, INFO )
      CALL SGEQRF( T2, T2, Q2, LDQ2, WORK, WORK( T2 + 1 ), 
     $             LWORK-T2, INFO )

      S_QR2 = SECOND() - S_QR2
*
*     8.4. Make sure the diagonal entries are positive, 
*          Save additional sign info in WORK( q_2 - R+1, 2*q_2 - R )
*      
        DO 190 I = 1, T2
          IF( Q2( I, I ) .LT. ZERO ) THEN
            CALL SSCAL( T2-I+1, -ONE, Q2( I, I ), LDQ2 )
            WORK( T2+I ) = -ONE
          ELSE
            WORK( T2+I ) = ONE
          END IF
  190   CONTINUE
*
*     8.5. Compute V = V * Ur2, 
*          where Ur2 is the Houlseholder reflectors.
*
*        T = MIN( 1+R - L + q_2 , 1+R )
        T = 1 + R - MAX( L - q_2, 0 )          
*        
        CALL SORM2R( 'R', 'N', P, T2, T2, Q2, LDQ2, WORK, V( 1, T ),
     $               LDV, WORK( 2*T2+1 ), INFO )        
*
*     8.6. Change the signs of V(: q_1+1-T:q_1 ) according to 
*          WORK( T+1: 2T )
*
      DO 210 J = 1, T2
        IF( WORK( T2 + J ) .LT. ZERO ) THEN
          DO 200 I = 1, P 
            V( I, T+J-1 ) = -V( I, T+J-1 )
  200     CONTINUE
        END IF
  210 CONTINUE

      END IF
*  
*=============== PRINT RESULT AFTER STEP 8 =================
*
*      print *, " % after STEP 8:"
*      print *, " V = ..."      
*      CALL PTMATRIX( P, P, V, LDV )
*===========================================================
*
*     
  900 CONTINUE
       S_ALL = SECOND() - S_ALL
*      PRINT *, "in SORCS2.f, tIME: "
*      WRITE( *, FMT = 9990 ) M, P, L, S_SVD1, S_SVD2, 
*     $     (100*( S_SVD1 + S_SVD2 ) / S_ALL),  S_QR1, S_QR2, 
*     $     (100*(S_QR1 + S_QR2)/S_ALL), S_ALL

 9990  FORMAT(  I4, ' & ', I4, ' & ', I4, 
     $          ' & ',  F8.5, ' & ',  F8.5, 
     $          ' & ',  F8.2, ' & ',  F8.5, 
     $          ' & ',  F8.5, ' & ',  F8.2, ' & ', F8.5, ' \\') 
*
      RETURN
*
*     End of SORCS2 
*
      END      
