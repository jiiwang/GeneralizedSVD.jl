      SUBROUTINE SCKCSD( NM, MVAL, PVAL, NVAL, ISEED, THRESH,
     $                   NMAX, AB, A, AF, B, BF, U, V, Q, ALPHA, BETA, 
     $                   WORK, RWORK, NIN, NOUT, INFO )
*
*  -- LAPACK test routine (version 0.1) --
*     Univ. of California Davis
*     July, 2005
*      
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            INFO, NIN, NM, NOUT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 ), MVAL( * ), NVAL( * ), NMAX,
     $                   PVAL( * )
      REAL               A( * ), AF( * ), B( * ), BF( * ), 
     $                   ALPHA( * ), BETA( * ), AB( 2*NMAX, * ),
     $                   Q( * ), U( * ), V( * ), WORK( * ), RWORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SCKCSD tests SORCSD:
*         the CSD for M-by-N matrix Q1 and P-by-N matrix Q2.
*
*  Arguments
*  =========
*
*  NM      (input) INTEGER
*          The number of values of M contained in the vector MVAL.
*
*  MVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix row dimension M.
*
*  PVAL    (input) INTEGER array, dimension (NP)
*          The values of the matrix row dimension P.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix column dimension N.
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator.  The array
*          elements should be between 0 and 4095, otherwise they will be
*          reduced mod 4096, and ISEED(4) must be odd.
*          On exit, the next seed in the random number sequence after
*          all the test matrices have been generated.
*
*  THRESH  (input) REAL
*          The threshold value for the test ratios.  A result is
*          included in the output file if RESULTS >= THRESH.  To have
*          every test ratio printed, use THRESH = 0.
*
*  AB      (workspace) REAL array, dimension (2NMAX, NMAX)
*
*  A       (workspace) REAL array, dimension (NMAX*NMAX)
*
*  AF      (workspace) REAL array, dimension (NMAX*NMAX)
*
*  B       (workspace) REAL array, dimension (NMAX*NMAX)
*
*  BF      (workspace) REAL array, dimension (NMAX*NMAX)
*
*  U       (workspace) REAL array, dimension (NMAX*NMAX)
*
*  V       (workspace) REAL array, dimension (NMAX*NMAX)
*
*  Q       (workspace) REAL array, dimension (NMAX*NMAX)
*
*  ALPHA   (workspace) REAL array, dimension (NMAX)
*
*  BETA    (workspace) REAL array, dimension (NMAX*NMAX)
*
*  WORK    (workspace) REAL array, dimension (NMAX*NMAX)
*
*  RWORK   (workspace) REAL array, dimension (NMAX*NMAX)
*
*  NIN     (input) INTEGER
*          The unit number for input.
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  INFO    (output) INTEGER
*          = 0 :  successful exit
*          > 0 :  If SLATMS returns an error code, the absolute value
*                 of it is returned.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 5 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRSTT
      CHARACTER          DISTA, DISTB
      CHARACTER*3        PATH
      INTEGER            I, IINFO, IM, LDA, LDAB, LDWZ,
     $                   LDB, LDQ, LDR, LDU, LDV, LWORK, M, 
     $                   N, NFAIL, NRUN, NT, P
*     ..
*     .. Local Arrays ..
      REAL               RESULTS( NTESTS )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAHDG, ALASUM, SLAROR, SCSDTS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      print *, "Used SCKCSD ready on July 7, 05"
      PATH( 1: 3 ) = 'CSD'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = .TRUE.
      LDA = NMAX
      LDB = NMAX
      LDAB = 2*NMAX
      LDU = NMAX
      LDV = NMAX
      LDQ = NMAX
      LDR = NMAX
      LWORK = NMAX*NMAX
      LDWZ = NMAX
            print *, "PATH = ", PATH, 
     $               " NIN, NOUT ",      NIN, NOUT      
*
*     Do for each value of M in MVAL.
*
      DO 30 IM = 1, NM
         M = MVAL( IM )
         P = PVAL( IM )
         N = NVAL( IM )
*
*        Generate M+P by N matrix AB
*
         CALL SLAROR( 'Lift or Right', 'Init', M+P, N, AB, LDAB, 
     $          ISEED, WORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            WRITE( NOUT, FMT = 9999 )IINFO
            INFO = ABS( IINFO )
            GO TO 30
         END IF
         CALL SLACPY( 'All', M, N, AB, LDAB, A, LDA )
         CALL SLACPY( 'All', P, N, AB( M + 1, 1 ), LDAB, B, LDB )
*
         CALL SCSDTS( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU,
     $                V, LDV, Q, LDQ,
     $                ALPHA, BETA, WORK, LWORK, RWORK, RESULTS )
*
*        Print information about the tests that did not
*        pass the threshold.
*
         DO 10 I = 1, NTESTS
            IF( RESULTS( I ).GE.THRESH ) THEN
               IF( NFAIL.EQ.0 .AND. FIRSTT ) THEN
                  FIRSTT = .FALSE.
                  CALL ALAHDG( NOUT, PATH )
               END IF
               WRITE( NOUT, FMT = 9998 )M, P, N, 'Orthonormal', I,
     $            RESULTS( I )
               NFAIL = NFAIL + 1
            END IF
   10    CONTINUE
         NRUN = NRUN + NTESTS
   30 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, 0 )
*
 9999 FORMAT( ' SLATMS in SCKGSV   INFO = ', I5 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', N=', I4, ', type ', A,
     $      ', test ', I2, ', ratio=', G13.6 )
      RETURN
*
*     End of SCKCSD
*
      END
