      SUBROUTINE SGGPSV( JOBU, JOBVT, M, K, N, A, LDA, B, LDB,
     $                   S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
      IMPLICIT  NONE
*
*  -- LAPACK driver routine ( version 0.1 ) --
*     Univ. of California Davis,
*     April, 2005
*      
*     .. Scalar Arguments ..
      CHARACTER          JOBU, JOBVT
      INTEGER            M, K, N, LDA, LDB, LDU, LDVT, INFO, LWORK
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ),
     $                   S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SGGPSV computes the produce singular value decomposition (PSVD)
*  of an M-by-N real matrix A and N-by-P real matrix B:
*
*      A * B = U * SIGMA * V'
*
*  where SIGMA is the diagonal matrix that contains the singular 
*  values of the product A*B, U and V are orthogonal matrices.
*  Z' denotes the transpose of Z. 
*  The first min(M, N) columns of U and V are the left and right
*  singular vectors of A*B.
*
*  Note that the routine returns V', not V.      
*  
*  Arguments
*  =========
*
*  JOBU    (input) CHARACTER*1
*          = 'U':  Orthogonal matrix U is computed;
*          = 'N':  U is not computed.
*
*  JOBVT    (input) CHARACTER*1
*          = 'V':  Orthogonal matrix VT is computed;
*          = 'N':  VT is not computed.
*
*  M       (input) INTEGER
*          The number of rows of the input matrix A.  M >= 0.
*
*  K       (input) INTEGER
*          The number of columns of the input matrix A.  
*          The number of rows of the input matrix B.  K >= 0.      
*
*  N       (input) INTEGER
*          The number of columns of the matrix B. N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,K)
*          On entry, the M-by-K matrix A.
*          On exit, entries are those as returned from SGGBRD.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input/output) REAL array, dimension (LDB,N)
*          On entry, the K-by-N matrix B. 
*          On exit, entries are those as returned from SGGBRD.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,K).
*
*  S       (output) REAL array, dimension (min(M, N))
*          The singular values of A*B, sorted so that S(i) >= S(i+1).
*
*  U       (output) REAL array, dimension (LDU, N)
*          If JOBU = 'U', U contains the M-by-M orthogonal matrix U.
*          If JOBU = 'N', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U.
*          If JOBU = 'U', LDU >= M.
*          Otherwise, LDU >= 1
*
*  VT      (output) REAL array, dimension (LDVT, N)
*          If JOBV = 'VT', VT contains the N-by-N orthogonal matrix VT.
*          If JOBV = 'N',  VT is not referenced.
*
*  LDVT    (input) INTEGER
*          The leading dimension of the array VT. 
*          If JOBVT = 'V', LDVT >= N.
*          Otherwise, LDVT >= 1.
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*
*  LWORK   (input) INTEGER 
*          The dimension of the array WORK.  LWORK >= 1
*          LWORK >= max( M, N, K ) + 4*min( M, N )
*          Note: LWORK is somewhat conservative.
*      
*  INFO    (output)INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )      
*     ...      
*     .. Local Scalars ..
      LOGICAL            WANTU, WANTVT 
      INTEGER            ITAUP, ITAUQ, IE, IP, IQ, 
     $                   IMIN, WK, IWK
      REAL               S_ALL, S_BRD, S_SQR, S_GEN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SECOND
      EXTERNAL           LSAME, SECOND
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGGBRD, SLASET, SLACPY, SORGQR, SORGLQ, 
     $                   SBDSQR, XERBLA 
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     Executable Statements ..
*
*     Test the input parameters
*
      WANTU  = LSAME( JOBU,  'U'  )
      WANTVT = LSAME( JOBVT, 'V' )
*
      INFO = 0
      IF( .NOT.( WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
        INFO = -1
      ELSE IF( .NOT.( WANTVT .OR. LSAME( JOBVT, 'N' ) ) ) THEN
        INFO = -2
      ELSE IF( M.LT.0 ) THEN
        INFO = -3
      ELSE IF( K.LT.0 ) THEN
        INFO = -4
      ELSE IF( N.LT.0 ) THEN
        INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
        INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
        INFO = -9
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
        INFO = -12
      ELSE IF( LDVT.LT.1 .OR. ( WANTVT .AND. LDVT.LT.N ) ) THEN
        INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'SGGPSV', -INFO )
        RETURN
      END IF
*
*     Initialize integer variables
*   
*
      S_ALL = SECOND( )
*      
      IMIN  = MIN( M, N )
      IQ = MIN( M-1, K, N   )
*      IP = MIN( M  , K, N-2 )
      IP = MIN( M  , K, N-1 )      
*
      IE = 1
      ITAUP = IE + IMIN
      ITAUQ = ITAUP + IMIN
      IWK   = ITAUQ + IMIN      
      WK = LWORK - 3*IMIN - 1
*     
*     Bidiagonalize A*B without explicitly multipling them out.
*     Place the diagonal entries in S and and off-diagonal in
*     WORK( IE:IMIN-1 )
*
      S_BRD = SECOND( )       
*    
      CALL SGGBRD( M, K, N, A, LDA, B, LDB, S, WORK( IE ),
     $             WORK( ITAUQ ), WORK ( ITAUP ),  WORK( IWK ),
     $             WK, INFO )
*
       S_BRD = SECOND( ) - S_BRD       
*      
      IF( WANTU .AND. WANTVT ) THEN
*
*     Place the reflectors from A to U, and from B to VT
*      
       CALL SLASET( 'All', M, M, ZERO, ZERO, U, LDU )
       CALL SLASET( 'All', N, N, ZERO, ZERO, VT, LDVT )
       CALL SLACPY( 'Lower', M-1, IQ, A( 2, 1 ), LDA, U( 2, 1 ), LDU )
       CALL SLACPY( 'Upper', IP, N-1, B( 1, 2 ), LDB, VT( 1, 2 ), LDVT)
*
*     
        S_GEN = SECOND() 
*
*       Use SORGBR to generate orthogonal matrices, use work space
*       min(M,N)
*
        CALL SORGBR( 'Q', M, M, N, U, LDU, WORK( ITAUQ ), WORK( IWK ),
     $               WK, INFO )
        CALL SORGBR( 'P', N, N, M, VT, LDVT, WORK( ITAUP ), WORK( IWK ),
     $               WK, INFO )
* 
        S_GEN = SECOND( ) - S_GEN 
*          
        IQ = M
        IP = N
*
      ELSE IF( WANTU ) THEN
*
*       Place the reflectors from A to U only
*      
        CALL SLASET( 'All', M, M, ZERO, ZERO, U, LDU )
        CALL SLACPY( 'Lower', M-1, IQ, A( 2, 1 ), LDA, U( 2, 1 ), LDU )
*
*       Generate the orthogonal matrix U only
*      
        CALL SORGBR( 'Q', M, M, N, U, LDU, WORK( ITAUQ ), WORK( IWK ),
     $               WK, INFO )
        IQ = M
        IP = 0
*        
      ELSE IF( WANTVT ) THEN
*
*       Place the reflectors from B to VT only
*      
        CALL SLASET( 'All', N, N, ZERO, ZERO, VT, LDVT )
        CALL SLACPY( 'Upper', IP, N-1, B( 1, 2 ), LDB, VT( 1, 2 ), LDVT)
*
*       Generate the orthogonal matrix VT only
*      
        CALL SORGBR( 'P', N, N, M, VT, LDVT, WORK( ITAUP ), WORK( IWK ),
     $               WK, INFO )
        IQ = 0
        IP = N
*        
      ELSE
*
*       Both U and VT are not computed
*        
        IQ = 0
        IP = 0
*        
      END IF
*
*     Only WORK(1:IMIN-1) is occupied with sub-diagional entries.
*      
      IWK = IMIN
*           
      S_SQR = SECOND()
*
*     USE SBDSQR to compute the singular alues and vectors, 
*     use 4*min(M,N) work space
*      
      IF( M .LT. N ) THEN
*
*       M >= N, compute the SVD of a lower bidiagonal matrix
*
        CALL SBDSQR( 'Lower', IMIN, IP, IQ, 0, S, WORK( IE ), 
     $                VT, LDVT, U, LDU, A, LDA, WORK ( IWK ), INFO )
      ELSE
*
*       M < N, compute the SVD of an upper bidiagonal matrix
*         
        CALL SBDSQR( 'Upper', IMIN, N, M, 0, S, WORK( IE ),
     $                VT, LDVT, U, LDU, A, LDA, WORK ( IWK ), INFO ) 
         
      END IF        
*
      S_SQR = SECOND( ) - S_SQR
      S_ALL = SECOND( ) - S_ALL
      print *, "In SGGPSV, time: "
      WRITE( *, FMT=9999 ) M, K, N, 
     $  S_BRD, 100*S_BRD/S_ALL, S_GEN, 100*S_GEN/S_ALL, 
     $  S_SQR, 100*S_SQR/S_ALL, S_ALL
*
 9999 FORMAT( I4, ' & ', I4, ' & ', I4, ' & ', 
     $        F8.5, ' & ', F8.2, ' & ', 
     $        F8.5, ' & ', F8.2,
     $        F8.5, ' & ', F8.2, ' & ', F8.5, ' \\ ')
*      
      RETURN
*
*     End of SGGPSV
*
      END
