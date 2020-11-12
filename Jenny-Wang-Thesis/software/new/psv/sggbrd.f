      SUBROUTINE SGGBRD( M, K, N, A, LDA, B, LDB, D, E, TAUQ, TAUP,
     $                   WORK, LWORK,  INFO )
      IMPLICIT NONE
*
*  -- LAPACK driver routine ( version 0.1 ) --
*     Univ. of California Davis,
*     April, 2005
*      
*     .. Scalar Arguments ..
      INTEGER            M, N, K, LDA, LDB, LWORK, INFO
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ),
     $                   D( * ), E( * ), TAUP( * ), TAUQ( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SGGBRD reduces the product of two general matrices, A of size M-by-N
*  and B of size N-by-P, to upper bidiagonal form G by a sequence of 
*  orthogonal transformations:
*
*               Q' * (A * B) * P = G
*
*  where Q and P are orthogonal matrices, Z' denotes the transpose of Z.
*  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
*  
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows in the matrix A.  M >= 0.
*
*  K       (input) INTEGER
*          The number of columns in the matrix A. 
*          The number of rows in the matrix B.  K >= 0.
*       
*  N       (input) INTEGER
*          The number of columns in the matrix B.  N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA, K)
*          On entry, the M-by-K general matrix.
*          On exit,
*          if M >= N, the elements below the diagonal, 
*            with the array TAUQ, represent the orthogonal matrix Q 
*            as a product of elementary reflectors.
*          if M < N, the elements below the first subdiagonal,  
*            with the array TAUQ, represent the orthogonal matrix Q 
*            as a product of elementary reflectors. 
*          The rest of the entries in A are set to zeros.     
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (input/output) REAL array, dimension (LDB, K)
*          On entry, the K-by-N general matrix.
*          On exit,
*          if M >= N, the elements above the first superdiagonal,
*            with the array TAUP, represent the orthogonal matrix P as
*            a product of elementary reflectors.
*          if M < N, the elements above the diagonal 
*            with the array TAUP, represent the orthogonal matrix P as  
*            a product of elementary reflectors.      
*          The rest of the entries in B are set to zeros.
*          See Further Details.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  D       (output) REAL array, dimension (min(M,N))
*          The diagonal elements of the bidiagonal matrix G:
*
*  E       (output) REAL array, dimension (min(M,N) - 1)
*          if M >= N, the upper-diagonal elements of the bidiagonal 
*              matrix G,
*          If M < N, the lower-diagonal elements of the bidiagonal 
*              matrix G.      
*
*  TAUQ    (output) REAL array dimension (min(M,N))
*          The scalar factors of the elementary reflectors which
*          represent the orthogonal matrix Q. See Further Details.
*
*  TAUP    (output) REAL array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors which
*          represent the orthogonal matrix P. See Further Details.
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  
*          LWORK >= max( M, N, K ) + min( M, N ). 
*
*  INFO    (output) INTEGER
*          = 0:  successful exit 
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrices Q and P are represented as products of elementary
*  reflectors:
*
*  If M >= N,
*
*     Q = H(1) H(2) . . . H(t)  and  P = W(1) W(2) . . . W(r)
*     where t = min(m-1, k, n) and r = min(m, k, n-2).     
*
*  Each H(i) and W(i) has the form:
*
*     H(i) = I - tauq * v * v'  and W(i) = I - taup * u * u'
*
*  where tauq and taup are real scalars, and v and u are real vectors;
*  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
*  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in B(i,i+2:n);
*  tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  If M < N,
*     Q = H(1) H(2) . . . H(t)  and  P = G(1) G(2) . . . G(r)
*     where t = min(m, k-1) and r = min(m-2, k, n).
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are real scalars, and v and u are real vectors;
*  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
*  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in B(i,i+1:n);
*  tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  The contents of A and B on exit are illustrated by the following
*  examples:
*
*  M = 6 K = 7 and N = 5 ( M >= N )
*   A =                        B =  
*    ( -  -  -  -  -  -  - )       (  -  -  u1 u1 u1 )    
*    ( v1 -  -  -  -  -  - )       (  -  -  -  u2 u2 )    
*    ( v1 v2 -  -  -  -  - )       (  -  -  -  -  u3 )    
*    ( v1 v2 v3 -  -  -  - )       (  -  -  -  -  -  )    
*    ( v1 v2 v3 v4 -  -  - )       (  -  -  -  -  -  )
*    ( v1 v2 v3 v4 v5 -  - )       (  -  -  -  -  -  )
*                                  (  -  -  -  -  -  )
*
*  M = 5, K = 7, N = 6   ( M < N )
*   A =                        B =  
*    ( -  -  -  -  -  -  - )       (  -  u1  u1  u1  u1  u1 )    
*    ( -  -  -  -  -  -  - )       (  -  -   u2  u2  u2  u2 )    
*    ( v1 -  -  -  -  -  - )       (  -  -   -   u3  u3  u3 )    
*    ( v1 v2 -  -  -  -  - )       (  -  -   -   -   u4  u4 )    
*    ( v1 v2 v3 -  -  -  - )       (  -  -   -   -   -   u5 )
*                                  (  -  -   -   -   -   -  )
*                                  (  -  -   -   -   -   -  )
*    
*
*  where vi denotes an element of the vector defining H(i), 
*  and ui an element of the vector defining W(i), 
*  and the dash denotes an zero element.     
*  The main diagonal entries of A*B are stored in D and the 
*  off-diagonal entries of A*B are stored in E.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )      
*     ...      
*     .. Local Scalars ..
      INTEGER            I, J, Ia, Ib, Iab, Imin
      REAL               TMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARF, SLARFG, SLASET, SGEMV, SCOPY, SSDOT, 
     $                   XERBLA

*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( K.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, K ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGGBRD', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*      
      J = MAX( MAX( M, K ), N )
      IF( J .EQ. 0 ) THEN
         RETURN
      ELSE IF ( J .EQ. 1 ) THEN
         IF ((M .EQ. 1) .AND. (K .EQ. 1) .AND. (N .EQ. 1 )) THEN
           D( 1 ) = A( 1, 1 ) * B( 1, 1 )
           A( 1, 1 ) = ONE
           B( 1, 1 ) = ONE
           RETURN
         ELSE 
           RETURN
         END IF
      END IF
*      
*     Initialize output variables
*      
      Imin = MIN( M, N )
      CALL SLASET('All', Imin,   1, ZERO, ZERO, TAUQ, 1 )
      CALL SLASET('All', Imin,   1, ZERO, ZERO, TAUP, 1 )
      CALL SLASET('All', Imin,   1, ZERO, ZERO,    D, 1 )
      CALL SLASET('All', Imin-1, 1, ZERO, ZERO,    E, 1 )
*
      IF( M .LT. N ) THEN
        GOTO 100
      END IF
*      print *, "%in M >= N"
*
*     This is the case where M >= N
*      
      Ia  = MIN( M-1, MIN( K,   N ) )
      Ib  = MIN( M-1, MIN( K-1, N ) )
      Iab = MIN( M  , MIN( K, N-2 ) )
*      
      J = MAX( Ia, Iab )
      DO I = 1, J
         IF( I .LE. Ib ) THEN
*
*          Form the Householder reflector H that annihilates B(I:K, I).
*
           CALL SLARFG( K-I+1, B( I, I ), B( I+1, I ), 1, TAUQ( I ))
*
*          Apply H' to the left of B
*
           TMP = B( I, I )
           B( I, I ) = ONE
           CALL SLARF( 'Left', K-I+1, N-I, B( I, I ), 1, TAUQ( I ), 
     $                 B( I, I+1 ), LDB, WORK )
*
*          Apply H to the right of A
*               
           CALL SLARF( 'Right', M, K-I+1, B( I, I ), 1, TAUQ( I ), 
     $                  A( 1, I ), LDA, WORK )
           B( I, I ) = TMP
*
*          Clear B( I+1:K, I )
*             
           CALL SLASET('All', K-I, 1, ZERO, ZERO, B( I+1, I ), 1 )
*         
         END IF
*         
         IF( I .LE. Ia ) THEN
*
*           Form the Householder reflector (H2) That annihilates A(I:M, I).
*
            CALL SLARFG( M-I+1, A( I, I ), A( I+1, I ), 1, TAUQ( I ))
*
*           Apply H2' to the left of A
*
            TMP = A( I, I )
            A( I, I ) = ONE
            CALL SLARF( 'Left', M-I+1, K-I, A( I, I ), 1, TAUQ( I ), 
     $                  A( I, I+1 ), LDA, WORK )
            A( I, I ) = TMP
*
         END IF

*
*        Multiply A(I, I:K) B(I:K, I:N), place result in WORK(1:N-I+1)
*
         CALL SGEMV( 'Transpose', K-I+1, N-I+1, ONE, B( I, I ), LDB, 
     $                  A( I, I ), LDA, ZERO, WORK, 1 )
*
         D( I ) = WORK( 1 )
         B( I, I+1 ) = ZERO
*
         IF( I .LE. Iab ) THEN
*
*           Form the Householder reflector (H0) that annihilates
*           WORK(2:N-I+1)
*
            CALL SLARFG( N-I, WORK( 2 ), WORK( 3 ), 1, TAUP( I ))
*
*           Apply H0 to the right of B(I+1:K, I+1:N)
* 
            E( I ) = WORK( 2 )
            WORK( 2 ) = ONE
            CALL SLARF( 'Right', K-I, N-I, WORK( 2 ), 1, TAUP( I ), 
     $                  B( I+1, I+1 ), LDB, WORK( N+1 ) )
*  
*           Save Householer reflector in B(I, I+2:N) 
*           
            CALL SCOPY( N-I-1, WORK( 3 ), 1, B( I, I+2 ), LDB )
*
          ELSE IF( I .LT. N ) THEN            
            E( I ) = WORK( 2 )
*
         END IF
*
      END DO
*            
      IF( I-1 .LT. Imin ) THEN
*
*       D( I ) = A( I, I:K ) * B( I:K, I )
*        
        CALL SSDOT( K-I+1, D( I ), A( I, I ), LDA, B( I, I ), 1 )
*        
      END IF
*
*     Clean up      
*      
      CALL SLASET('Upper', M, K, ZERO, ZERO, A(1, 1), LDA)
      CALL SLASET('Lower', K, N, ZERO, ZERO, B(1, 1), LDB)
      CALL SLASET('Lower',   M-I+1, K-I+1, ZERO, ZERO, A(I, I), LDA )
      CALL SLASET('Upper',   K-I+1, N-I+1, ZERO, ZERO, B(I, I), LDB )
*
*     Finish
*      
      GOTO 200
*
  100 CONTINUE
*      print *, "% in M < N "
*
*     This is the case where M < N
*      
      Ia  = MIN( M, K-1 )
      Ib  = MIN( M, K )
      Iab = MIN( M-2, K, N )
*
      J = MAX( Ib, Iab )
      DO I = 1, J
         IF( I .LE. Ia ) THEN
*
*          Form the Householder reflector H that annihilates A(I, I:K).
*
           CALL SLARFG( K-I+1, A( I, I ), A( I, I+1 ), LDA, TAUP( I ))
* 
*          Note: the vector v in the representation of H must
*                be copied to WORK first because it is a row vector
*                and the increment between elements of v in A is LDA
*                ( won't work correctly if use the v in A ).           
*
           WORK( 1 ) = ONE
           CALL SCOPY( K-I, A( I, I+1 ), LDA, WORK( 2 ), 1 )
*
*          Apply H to the right of A(I+1:M, I:K)
*
           CALL SLARF( 'Right', M-I, K-I+1, WORK, 1, TAUP( I ),
     $                 A( I+1, I ), LDA, WORK( K + 1 ) )
*
*          Apply H to the left of B(I:K, 1:N)
*
           CALL SLARF( 'Left', K-I+1, N, WORK, 1, TAUP( I ),
     $                 B( I, 1 ), LDB, WORK( K + 1 ) )         
*           A( I, I ) = TMP
*
*          Clear A( I, I+1:K )
*             
           CALL SLASET('All', 1, K-I, ZERO, ZERO, A( I, I+1 ), LDA )
*           
         END IF
*
         IF( I .LE. Ib ) THEN
*
*           Form the Householder reflector (H2) That annihilates B(I, I:N).
*
            CALL SLARFG( N-I+1, B( I, I ), B( I, I+1 ), LDB, TAUP( I ))
*
*           Again, copy v in the representation of H2 to WORK first
*            
           WORK( 1 ) = ONE
           CALL SCOPY( N-I, B( I, I+1 ), LDB, WORK( 2 ), 1 )
*
*           Apply H2 to the right of B(I+1:K, I:N)
*            
            CALL SLARF( 'Right', K-I, N-I+1, WORK, 1, TAUP( I ), 
     $                  B( I+1, I ), LDB, WORK( N+1 ) )
*           
         END IF         
*
*        Multiply A(I:M, I:K) B(I:K, I), place result in WORK( 1:M-I+1 )       
*
         CALL SGEMV( 'No Transpose', M-I+1, K-I+1, ONE, A( I, I ), LDA, 
     $                  B( I, I ), 1, ZERO, WORK, 1 )
*
         D( I ) = WORK( 1 )
         A( I+1, I ) = ZERO
         A( I, I ) = ZERO
*
         IF( I .LE. Iab ) THEN
*
*           Form the Householder reflector (H0) that annihilates WORK( 2:M-I+1 )
*
            CALL SLARFG( M-I, WORK( 2 ), WORK( 3 ), 1, TAUQ( I ))
*
*           Apply H0 to the left of A(I+1:M, I+1:K)
* 
            E( I ) = WORK( 2 )
            WORK( 2 ) = ONE
            CALL SLARF( 'Left', M-I, K-I, WORK( 2 ), 1, TAUQ( I ), 
     $                  A( I+1, I+1 ), LDA, WORK( M+1 ) )
*  
*           Save Householer reflector in A( I+2:M, I ) 
*           
            CALL SCOPY( M-I-1, WORK( 3 ), 1, A( I+2, I ), 1 )
*
          ELSE IF( I .LT. M ) THEN
            E( I ) = WORK( 2 )
          END IF
      END DO 
*
*     Clean up      
*      
*      CALL SLASET('Upper', M, K, ZERO, ZERO, A(1, 1), LDA)
      CALL SLASET('Lower', K, N, ZERO, ZERO, B(1, 1), LDB)
*      CALL SLASET('Lower',   M-I+1, K-I+1, ZERO, ZERO, A(I, I), LDA )
      CALL SLASET('Upper',   K-I+1, N-I+1, ZERO, ZERO, B(I, I), LDB )
       
*           
  200 RETURN
*
*     End of SGGBRD
*
      END
