      SUBROUTINE SPSVTS( M, P, N, A, AF, LDA, B, BF, LDB, C, 
     $                   U, LDU, Z, LDZ, ALPHA, WORK, 
     $                   LWORK, RESULTS )
        IMPLICIT NONE
*
*  -- LAPACK test routine (version 0.1) --
*     Univ of California Davis,
*     June, 2005
*     ..
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDU, LDZ, M, P, N, LWORK
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), AF( LDA, * ), C( LDA, * ), 
     $                   B( LDB, * ), BF( LDB, * ), 
     $                   U( LDU, * ), Z( LDZ, * ),
     $                   ALPHA( * ),  
     $                   RESULTS( 6 ), WORK( LWORK )        
*     ..
*
*  Purpose
*  =======
*
*  SPSVTS checks the product singular value decomposition (PSVD) 
*  routines.
*
*  SGGBRD reduces a pair of real general matrices, A (m-by-p) and 
*  B (p-by-n), to upper or lower bidiagonal form G by an orthogonal 
*  transformation:  Q' * A * B * P = G
*  (or A * B = Q * G * P').  The matrix G is upper bidiagonal if m >= n
*  and lower bidiagonal if m < n.   Z' denotes the transpose of Z.
*
*  SORGBR generates the orthogonal matrices Q and P' from SGGBRD.
*  Note that Q and P are not necessarily square.
*
*  SPSVTS computes the product singular value decomposition by first 
*  calling SGGBRD to reduce A*B to bidiagonal form followed by 
*  calling SBDSQR to find the singular value decomposition of G 
*  (U1 S V1' = G ).  The result of the product singular value 
*  decomposition is then obtained by combining the resulted orthogonal
*  matrices to form the final decomposed components: 
*     A * B = (Q U1) S (V1' P') = U S V' where 
*     U = Q U1 and V = P1 V1
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  P       (input) INTEGER
*          The number of columns of the matrix A.  P >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices B.  N >= 0.
*
*  A      (input) REAL array, dimension ( LDA, P )
*          The M-by-A matrix A
*
*  AF      (workspace) REAL array, dimension ( LDA, P )
*
*  LDA    (input) INTEGER
*          The leading dimension of the arrays A
*          LDA >= max( 1, M ).
*
*  B      (input) REAL array, dimension ( LDB, N )
*          On entry, the P-by-N matrix B.
*
*  BF      (workspace) REAL array, dimension ( LDB, N )
*
*  C      (workapce) REAL array, dimension ( LDA, N )
*          On exit, the M-by-N matrix C.
*
*  LDB    (input) INTEGER
*          The leading dimension of the arrays B and LDB
*          LDB >= max( 1, P ).
*
*  U       (output) REAL array, dimension( LDU, M )
*          The M by M orthogonal matrix U.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U. LDU >= max( 1, M ).
*
*  Z       (output) REAL array, dimension( LDV, N )
*          The N by N orthogonal matrix Z.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z. LDZ >= max( 1, N ).
*
*  ALPHA   (output) REAL array, dimension ( MIN( M, N ) )
*          The singular value pairs of the product A*B, the
*        
*  WORK    (workspace) REAL array, dimension ( LWORK )
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK,
*          LWORK = NMAX*NMAX
*
*  RESULTS  (output) REAL array, dimension >=  5 
*          The test ratios:
*  Test for SGGPSV:
*
*  (1)    | A B - U diag(S) VT | / ( |A B| max(M,N) ulp )
*
*  (2)    | I - U'U | / ( M ulp )
*
*  (3)    | I - VT VT' | / ( N ulp )
* 
*  Test for SGGBRD on bidiagonal matrix G:
*
*  (4)    | G - Q diag(S) PT | / ( |G| max(M,N) ulp )
*
*  (5)    | I - Q'Q | / ( M ulp )
*
*  (6)    | I - PT PT' | / ( N ulp )
*            
*
*  =====================================================================
*
        REAL ONE, ZERO
        PARAMETER( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalar .. 
        INTEGER I, MNMIN, ITAUP, ITAUQ, IWK, INFO, 
     $          IQ, IP
        REAL    ULP, ULFL, DUMMY
*     ..
*     .. External Subroutines ..
      External  SLASET, SGGPSV, SGEMM, SLACPY,  
     $          SBDT01, SORT01, SORGBR
*     .. Intrinsic Functions ..
      INTRINSIC MAX, MIN
*
      CALL SLACPY( 'All', M, P, A, LDA, AF, LDA )
      CALL SLACPY( 'All', P, N, B, LDB, BF, LDB )
*      
      MNMIN = MIN( M, N )

*
*  3. C = A*B
       CALL SGEMM('No Transpose', 'No Transpose', M, N, P, ONE, A, 
     $            LDA, B, LDB, ZERO, C, LDA )
*
*  4. Call PSVD routine ( for exact work size, see SGGPSV )
*

*      print *, "AF = ..."
*      call ptmatrix( M, P, AF, LDA )
*      print *, "BF = ..."
*      call ptmatrix( P, N, BF, LDB )
       
      CALL SGGPSV( 'U', 'V', M, P, N, AF, LDA, BF, LDB, ALPHA, U, LDU,
     $              Z, LDZ, WORK, LWORK, INFO )
      
*      print *, "Uf = ..."
*      call ptmatrix( M, M, U, LDU )
*      print *, "ZfT = ..."
*      call ptmatrix( N, N, Z, LDZ )
*      print *, "ALPHAf = ..."
*      call ptmatrix( 1, MIN(M, N), ALPHA, 1 )
      
      
*
*  5. Test for SGGPSV on A and B:
*      5.1. Test the residual of AB-USV' (Need WORK to be sized M+N)
*
      CALL SBDT01( M, N, 0, C, LDA, U, LDU, ALPHA, DUMMY, Z, LDZ, WORK,
     $              RESULTS( 1 ) )
*
*      5.2 Test the orthogonality of U
*
       IF( LWORK.LT.N*(N+1) ) THEN
         PRINT *, "WARN: LWORK < N(N+1) (SORT01 for Columns) !!"
       END IF
       CALL SORT01( 'Columns', M, M, U, LDU, WORK, LWORK, RESULTS( 2 ) )
*   
*      5.3 Test the orthogonality of Z'
*
       IF( LWORK.LT.M*(M+1)) THEN
         PRINT *, "WARN: LWORK < M*(M+1) (SORT01 for Rows) !!"
       END IF
*
       CALL SORT01( 'Rows', N, N, Z, LDZ, WORK, LWORK, RESULTS( 3 ) )
*
*=== Set up test for SGGBRD ===
*
*  4. Call SGGBRD routine ( for exact work size, see SGGBRD )
*
      I = 1
      ITAUP = I + MNMIN
      ITAUQ = ITAUP + MNMIN
      IWK   = ITAUQ + MNMIN
*     
       CALL SGGBRD( M, P, N, A, LDA, B, LDB, ALPHA, WORK( I ),
     $             WORK( ITAUQ ), WORK ( ITAUP ),  WORK( IWK ),
     $             LWORK-3*MNMIN-1, INFO )
*
*  5. Test for SGGBRD on G: 

*     5.1 Place Q and PT in U and VT respectively 
*
       IQ = MIN( M-1, P, N   )
       IP = MIN( M  , P, N-1 )
       CALL SLASET( 'All', M, M, ZERO, ZERO, U, LDU )
       CALL SLASET( 'All', N, N, ZERO, ZERO, Z, LDZ )
*       CALL SLACPY( 'Lower', M-1, IQ, A( 2, 1 ), LDA, U( 2, 1 ), LDU )
*       CALL SLACPY( 'Upper', IP, N-1, B( 1, 2 ), LDB, Z( 1, 2 ), LDZ )
       CALL SLACPY( 'Lower', M-1, IQ, A( 2, 1 ), LDA, U( 2, 1 ), LDU )
       CALL SLACPY( 'Upper', IP, N-1, B( 1, 2 ), LDB, Z( 1, 2 ), LDZ )

       CALL SORGBR( 'Q', M, M, N, U, LDU, WORK( ITAUQ ), WORK( IWK ),
     $               LWORK, INFO )
       CALL SORGBR( 'P', N, N, M, Z, LDZ, WORK( ITAUP ), WORK( IWK ),
     $               LWORK, INFO ) 
*
*     5.2 Test the residual of G-USV' (Need WORK to be sized M+N)
*
       CALL SBDT01( M, N, 1, C, LDA, U, LDU, ALPHA, WORK( 1 ), Z, LDZ, 
     $              WORK( IWK ), RESULTS( 4 ) )
*
*     5.3 Test the orthogonality of Q
*
       CALL SORT01( 'Columns', M, M, U, LDU, WORK, N*(N+1),RESULTS( 5 ))
*   
*     5.4 Test the orthogonality of P'
*
       CALL SORT01( 'Rows', N, N, Z, LDZ, WORK, M*(M+1), RESULTS( 6 ) )
*
      PRINT *, "in SPSVTS: stable"
      WRITE( *, FMT = 9990 ) M, P, N, RESULTS(1), 
     $     RESULTS(2), RESULTS(3), RESULTS(4), RESULTS(5), RESULTS(6)
       PRINT *, "------------------------------------------------"
*       pause
*
 9990  FORMAT( I4, ' & ', I4, ' & ', I4, ' & ', 
     $         F8.5, ' & ',  F8.5, ' & ',  F8.5, ' & ',  F8.5, 
     $         ' & ',  F8.5, ' & ',  F8.5, ' \\ ')        
*
       RETURN       
*       
*      End of SPSVTS
*
       END

