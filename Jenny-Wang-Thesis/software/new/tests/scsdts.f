      SUBROUTINE SCSDTS( M, P, K, Q1, F1, LDQ1, Q2, F2, LDQ2, U, LDU,
     $                   V, LDV, Z, LDZ, 
     $                   ALPHA, BETA, WORK, LWORK, RWORK, RESULTS )
*
*  -- LAPACK test routine (version 1.0) --
*     Univ. of California Davis,
*     June, 2005
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            LDQ1, LDQ2, LDU, LDV, LDZ, M, P, K, LWORK
*     ..
*     .. Array Arguments ..
      REAL               Q1( LDQ1, * ), F1( LDQ1, * ), Q2( LDQ2, * ),
     $                   F2( LDQ2, * ), U( LDU, * ), V( LDV, * ),
     $                   Z( LDZ, * ), ALPHA( * ), BETA( * ), 
     $                   RESULTS( 5 ), WORK( LWORK ), RWORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SCSDTS tests SORCSD, which computes the CSD of an M-by-N matrix Q1
*  and a P-by-N matrix Q2:
*               U'*Q1*Z = D1 and V'*Q2*Z = D2.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q1.  M >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix Q2.  P >= 0.
*
*  K       (input) INTEGER
*          The number of columns of the matrices Q1 and Q2.  K >= 0.
*
*  Q1      (input) REAL array, dimension ( LDQ1,K )
*          The M-by-K matrix Q1.
*
*  F1      (workspace) REAL array, dimension ( LDQ1, K )
*
*  LDQ1    (input) INTEGER
*          The leading dimension of the arrays Q1 and F1.
*          LDQ1 >= max( 1,M ).
*
*  Q2      (input) REAL array, dimension ( LDQ2, K )
*          On entry, the P-by-N matrix B.
*
*  F2      (workspace) REAL array, dimension ( LDQ2, K )
*
*  LDQ2    (input) INTEGER
*          The leading dimension of the arrays Q2 and F2.
*          LDQ2 >= max( 1,P ).
*
*  U       (output) REAL array, dimension( LDU, M )
*          The M by M orthogonal matrix U.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U. LDU >= max( 1, M ).
*
*  V       (output) REAL array, dimension( LDV,P )
*          The P by P orthogonal matrix V.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V. LDV >= max( 1, P ).
*
*  Z       (output) REAL array, dimension( LDZ, K )
*          The K by K orthogonal matrix Z.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z. LDZ >= max( 1, K ).
*
*  ALPHA   (output) REAL array, dimension ( K )
*  BETA    (output) REAL array, dimension ( K )
*          The cosine-side value pairs of Q1 and Q2, the
*          ``diagonal'' matrices D1 and D2 are constructed from
*          ALPHA and BETA, see subroutine SORCSD for details.
*
*  WORK    (workspace) REAL array, dimension ( LWORK )
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK,
*          LWORK >= MAX(3*MIN(N,K) + MAX(N, K), 5*MIN(N,K))
*                 + 2 MIN( M, P, K)
*
*  RWORK   (workspace) REAL array, dimension ( MAX( M, K, P ) )
*
*  RESULTS  (output) REAL array, dimension >=  5 
*          The test ratios:
*          RESULTS(1) = norm( U'*Q1*Z - D1 ) / ( MAX(M,K)*norm(Q1)*ULP)
*          RESULTS(2) = norm( V'*Q2*Z - D2 ) / ( MAX(P,K)*norm(Q1)*ULP)
*          RESULTS(3) = norm( I - U'*U ) / ( M*ULP )
*          RESULTS(4) = norm( I - V'*V ) / ( P*ULP )
*          RESULTS(5) = norm( I - Z'*Z ) / ( K*ULP )      
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, I, J
      REAL               Q1NORM, Q2NORM, ULP, UNFL, RESID
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY
      EXTERNAL           SLAMCH, SLANGE, SLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SORCSD, SLACPY, SLASET, SSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
      ULP = SLAMCH( 'Precision' )
      print *, "ULP = ", ULP
      UNFL = SLAMCH( 'Safe minimum' )
*
*     Copy the matrix Q1 to the array F1, Q2 to F2.
*
      CALL SLACPY( 'Full', M, K, Q1, LDQ1, F1, LDQ1 )
      CALL SLACPY( 'Full', P, K, Q2, LDQ2, F2, LDQ2 )
*
      Q1NORM = MAX( SLANGE( '1', M, K, Q1, LDQ1, WORK ), UNFL )
      Q2NORM = MAX( SLANGE( '1', P, K, Q2, LDQ2, WORK ), UNFL )
*
*     Factorize the matrices Q1 and Q2 in the arrays F1 and F2.
*      
      CALL SORCSD( 'Y', M, P, K, F1, LDQ1, F2, LDQ2,
     $             ALPHA, BETA, U, LDU, V, LDV, Z, LDZ, WORK, 
     $             LWORK, INFO)
      
*      print *, "Uf = ..."
*      call ptmatrix(M, M, U, LDU)
*      print *, "Vf = ..."
*      call ptmatrix(P, P, V, LDV)
*      print *, "Zf = ..."
*      call ptmatrix(K, K, Z, LDZ)
*      print *, "Alphaf = ..."
*      call ptmatrix(K, 1, ALPHA, 1)
*      print *, "Betaf = ..."
*      call ptmatrix(K, 1, BETA, 1)
*
*     Compute Q1:= U' Q1  ZT'
*
 
      CALL SGEMM( 'No Transpose', 'Transpose', M, K, K, ONE,
     $            Q1, LDQ1, Z, LDZ, ZERO, WORK, M )
*
      CALL SGEMM( 'Transpose', 'No transpose', M, K, M, ONE, U, LDU,
     $            WORK, M, ZERO, Q1, LDQ1 )
*
*     Compute Q1 := Q1 - D1
*
      
      DO 10 I = 1, MIN( M, K )
        Q1( I, I ) = Q1( I, I ) - ALPHA( I )
   10 CONTINUE
*
*     Compute norm( U'*Q1*Z - D1 ) / ( MAX(1,M,K)*norm(Q1)*ULP )
*     .
*
      RESID = SLANGE( '1', M, K, Q1, LDQ1, WORK )
*
      IF( Q1NORM.GT.ZERO ) THEN
        RESULTS( 1 ) = ( ( RESID / REAL( MAX(1,M,K) ) ) 
*     $                 ) / ULP
     $                 / Q1NORM ) / ULP
      ELSE
         RESULTS( 1 ) = ZERO
      END IF
*
*     Compute Q2 := V' Q2 ZT'
*
      CALL SGEMM( 'No transpose', 'Transpose', P, K, K, ONE, Q2, 
     $            LDQ2, Z, LDZ, ZERO, WORK, P )
*
      CALL SGEMM( 'Transpose', 'No transpose', P, K, P, ONE, V, LDV,
     $            WORK, P, ZERO, Q2, LDQ2 )
*
*     Compute Q2 := Q2 - D2
*
      J = MAX( K-P, 0 )
      DO 20 I = 1, MIN( K, P )
        Q2( I, I+J ) = Q2( I, I+J ) - BETA( I+J )
   20 CONTINUE
*
*     Compute norm( V'*Q2*Z - D2 ) / ( MAX(1,P,K)*norm(Q2)*ULP ) .
*
      RESID = SLANGE( '1', P, K, Q2, LDQ2, WORK )
      IF( Q2NORM.GT.ZERO ) THEN
        RESULTS( 2 ) = ( ( RESID / REAL( MAX(1,P,K ) ) )
*     $                 ) / ULP
     $                 / Q2NORM  ) / ULP
      ELSE
         RESULTS( 2 ) = ZERO
      END IF
*
*     Compute I - U'*U
*
      CALL SLASET( 'Full', M, M, ZERO, ONE, WORK, M )
      CALL SSYRK( 'Upper', 'Transpose', M, M, -ONE, U, LDU, ONE, WORK,
     $            M )
*
*     Compute norm( I - U'*U ) / ( M * ULP ) .
*
      RESID = SLANSY( '1', 'Upper', M, WORK, M, RWORK )
      RESULTS( 3 ) = ( RESID / REAL( MAX( 1, M ) ) ) / ULP
*
*     Compute I - V'*V
*
      CALL SLASET( 'Full', P, P, ZERO, ONE, WORK, P )
      CALL SSYRK( 'Upper', 'Transpose', P, P, -ONE, V, LDV, ONE, WORK,
     $            P )
*
*     Compute norm( I - V'*V ) / ( P * ULP ) .
*
      RESID = SLANSY( '1', 'Upper', P, WORK, P, RWORK )
      RESULTS( 4 ) = ( RESID / REAL( MAX( 1, P ) ) ) / ULP
*
*     Compute I - Z'*Z
*
      CALL SLASET( 'Full', K, K, ZERO, ONE, WORK, K )
      CALL SSYRK( 'Upper', 'No transpose', K, K, -ONE, Z, LDZ,
     $            ONE, WORK, K )
      RESID = SLANSY( '1', 'Upper', K, WORK, K, RWORK )
      RESULTS( 5 ) = ( RESID / REAL( MAX( 1, P ) ) ) / ULP
      CALL SORT01( 'Column', M, M, U, LDU, WORK, LWORK, RESULTS( 3 ) )
      CALL SORT01( 'Column', P, P, V, LDV, WORK, LWORK, RESULTS( 4 ) )
      CALL SORT01( 'Row', K, K, Z, LDZ, WORK, LWORK, RESULTS( 5 ) )
*
      PRINT *, "in SCSDTS, RESIDUALS: "
      WRITE( *, FMT = 9990 ) M, P, K, RESULTS(1), 
     $     RESULTS(2), RESULTS(3), RESULTS(4), RESULTS(5)
      
*       print *, "Q1NORM = ", Q1NORM, ", Q2NORM = ", Q2NORM
*       pause
*
 9990  FORMAT(  I4, ' & ', I4, ' & ', I4, ' & ', 
     $          ' & ',  F8.5, ' , & ',  F8.5, 
     $          ' & ',  F8.5, ' , & ',  F8.5, 
     $          ' & ',  F8.5, ' \\') 
*
*
*
      RETURN
*
*     End of SCSDTS
*
      END
