      SUBROUTINE SQSVTS( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU,
     $                   V, LDV, Q, LDQ, 
     $                   ALPHA, BETA, R, LDR,
     $                   IWORK, WORK, LWORK, RWORK,RESULTS )
      IMPLICIT NONE
*
*  -- LAPACK test routine (version 0.1) --
*     Univ. of California Davis,
*     June, 2005
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDU, LDV, LDQ, LDR, LWORK, M, P, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               A( LDA, * ), AF( LDA, * ), B( LDB, * ),
     $                   BF( LDB, * ), U( LDU, * ), V( LDV, * ),
     $                   Q( LDQ, * ), R( LDR, * ), ALPHA( * ),
     $                   BETA( * ), RESULTS( 5 ), RWORK( * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  SQSVTS tests SGGQSV, which computes the GSVD of an M-by-N matrix A
*  and a P-by-N matrix B:
*               U'*A*Q = D1*R and V'*B*Q = D2*R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B.  N >= 0.
*
*  A       (input) REAL array, dimension (LDA,M)
*          The M-by-N matrix A.
*
*  AF      (output) REAL array, dimension (LDA,N)
*          Details of the GSVD of A and B, as returned by SGGSVD,
*          see SGGSVD for further details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A and AF.
*          LDA >= max( 1,M ).
*
*  B       (input) REAL array, dimension (LDB,P)
*          On entry, the P-by-N matrix B.
*
*  BF      (output) REAL array, dimension (LDB,N)
*          Details of the GSVD of A and B, as returned by SGGSVD,
*          see SGGSVD for further details.
*
*  LDB     (input) INTEGER
*          The leading dimension of the arrays B and BF.
*          LDB >= max(1,P).
*
*  U       (output) REAL array, dimension(LDU,M)
*          The M by M orthogonal matrix U.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U. LDU >= max(1,M).
*
*  V       (output) REAL array, dimension(LDV,M)
*          The P by P orthogonal matrix V.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V. LDV >= max(1,P).
*
*  Q       (output) REAL array, dimension(LDQ,N)
*          The N by N orthogonal matrix Q.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= max(1,N).
*
*  ALPHA   (output) REAL array, dimension (N)
*  BETA    (output) REAL array, dimension (N)
*          The generalized singular value pairs of A and B, the
*          ``diagonal'' matrices D1 and D2 are constructed from
*          ALPHA and BETA, see subroutine SGGSVD for details.
*
*  R       (output) REAL array, dimension(LDQ,N)
*          The upper triangular matrix R.
*
*  LDR     (input) INTEGER
*          The leading dimension of the array R. LDR >= max(1,N).
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK,
*          LWORK >= 6*max(M,P,N)*max(M,P,N)+7 * max(M, 3N, P).
*
*  RWORK   (workspace) REAL array, dimension (max(M,P,N))
*
*  RESULTS  (output) REAL array, dimension (5)
*          The test ratios:
*          RESULTS(1) = norm( U'*A*Q - D1*R ) / ( MAX(M,N)*norm(A)*ULP)
*          RESULTS(2) = norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP)
*          RESULTS(3) = norm( I - U'*U ) / ( M*ULP )
*          RESULTS(4) = norm( I - V'*V ) / ( P*ULP )
*          RESULTS(5) = norm( I - Q'*Q ) / ( N*ULP )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, K, L, I, J, CT
      REAL               ANORM, BNORM, ULP, UNFL, RESID
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY
      EXTERNAL           SLAMCH, SLANGE, SLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SGGSVJ, SLACPY, SLASET, SSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      ULP = SLAMCH( 'Precision' )
      UNFL = SLAMCH( 'Safe minimum' )  
*
*     Copy the matrix A to the array AF.
*
      CALL SLACPY( 'Full', M, N, A, LDA, AF, LDA )
      CALL SLACPY( 'Full', P, N, B, LDB, BF, LDB )
*
      ANORM = MAX( SLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( SLANGE( '1', P, N, B, LDB, RWORK ), UNFL )
      CALL SLASET( 'ALL', N, N, ZERO, ZERO, R, LDR  )
*
*     Factorize the matrices A and B in the arrays AF and BF.
*
      CALL SGGQSV( 'U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB,
     $             ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, 
     $             IWORK, WORK, LWORK, INFO )
*
*     Copy R
*
      DO 20 I = 1, MIN( K+L, M )
         DO 10 J = I, K + L
            R( I, J ) = AF( I, N-K-L+J )
   10    CONTINUE
   20 CONTINUE
*
      IF( M-K-L.LT.0 ) THEN
*     Get B(1:K+L-M, N+M-K-L+1:N) to R(M+1:K+L, N+M-K-L+1:N)
         DO 40 I = M + 1, K + L
            DO 30 J = I, K + L
               R( I, J ) = BF( I-M, N-K-L+J )
   30       CONTINUE
   40    CONTINUE
      END IF
*
*      print *, "Rf = ..."
*      call ptmatrix( K+L, K+L, R, LDR )
*      print *, "Uf = ..."
*      call ptmatrix( M, M, U, LDU )
*      print *, "Vf = ..."
*      call ptmatrix( P, P, V, LDV )
*      print *, "Qf = ..."
*      call ptmatrix( N, N, Q, LDQ )
*      print *, "Alphaf = ..."
*      call ptmatrix( 1, N, ALPHA, 1 )
*      print *, "Betaf = ..."
*      call ptmatrix( 1, N, BETA, 1 )
*      print *, "K = ", K, "; L = ", L
*
*
*     Compute A:= U'*A*Q
*
      CALL SGEMM( 'No transpose', 'No transpose', M, N, N, ONE,
     $            A, LDA, Q, LDQ, ZERO, WORK, M )
*
      CALL SGEMM( 'Transpose', 'No transpose', M, N, M, ONE, U, LDU,
     $            WORK, M, ZERO, A, LDA )
*
*     Compute A := A - D1 R
*
      IF( P.GE.K+L ) THEN
          CT = K
      ELSE
         CT = K + L - P
*        IF( M.GE.K+L ) THEN
*          CT = K+L-P
*        ELSE
*          CT = K+L-P
*        END IF
      END IF

      DO 60 I = 1, CT
        DO 50 J = I, K + L
          A( I, N-K-L+J ) = A( I, N-K-L+J ) - R( I, J )
   50   CONTINUE
   60 CONTINUE
*     
      DO 80 I = CT+1, MIN( K+L, M )
        DO 70 J = I, K + L
          A( I, N-K-L+J ) = A( I, N-K-L+J ) - ALPHA( I )*R( I, J )
   70   CONTINUE
   80 CONTINUE 
*
      RESID = SLANGE( '1', K, N, A, LDA, RWORK )          
      RESID = SLANGE( '1', M-K, N, A, LDA, RWORK )
*
      IF( ANORM.GT.ZERO ) THEN
         RESULTS( 1 ) = ( ( RESID / REAL( MAX(1,M,N) ) ) / ANORM ) / ULP
      ELSE
         RESULTS( 1 ) = ZERO
      END IF
*
*     Compute B:= V'*B*Q      
*
      CALL SGEMM( 'No transpose', 'No transpose', P, N, N, ONE, B, LDB,
     $            Q, LDQ, ZERO, WORK, P )
*
      CALL SGEMM( 'Transpose', 'No transpose', P, N, P, ONE, V, LDV,
     $            WORK, P, ZERO, B, LDB )
*
*     Compute B := B - D2 * R
*
      DO I = 1, L
         DO J = 1, L 
            B( I, N - L + J ) = B( I, N - L + J ) - 
     $                          BETA( I + K ) * R( K + I, K + J )
         END DO
      END DO
*
*     Compute norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP ) .
*
      RESID = SLANGE( '1', P, N, B, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULTS( 2 ) = ( ( RESID / REAL( MAX(1,P,N ) ) )/BNORM ) / ULP
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
*     Compute I - Q'*Q
*
      CALL SLASET( 'Full', N, N, ZERO, ONE, WORK, N )
      CALL SSYRK( 'Upper', 'Transpose', N, N, -ONE, Q, LDQ, ONE, WORK,
     $            N )
*
*     Compute norm( I - Q'*Q ) / ( N * ULP ) .
*
      RESID = SLANSY( '1', 'Upper', N, WORK, N, RWORK )
      RESULTS( 5 ) = ( RESID / REAL( MAX( 1, N ) ) ) / ULP
*
      print *, "in SQSVTS, residuals: "
      WRITE( *, FMT=9998)M,P,N,RESULTS(1), RESULTS(2), RESULTS(3), 
     $     RESULTS(4), RESULTS(5)

 9998 FORMAT( I4, ' & ', I4, ' & ', I4, ' & ',F8.5, ' & ', F8.5, ' & ',
     $       F8.5, ' & ', F8.5, ' & ', F8.5, ' \\' )      
      RETURN
*
*     End of SQSVTS
*
      END
