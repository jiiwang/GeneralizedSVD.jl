      PROGRAM SCHKEE2
*
      IMPLICIT NONE
*
*  -- LAPACK test routine (version 0.1) --
*     Univ. of California Davis,
*     July, 2005
*
*  Purpose
*  =======
*
*  SCHKEE2 tests the REAL LAPACK subroutines CSD, QSV, PSV.
*
*  CSD (Cosin-Sine Decomposition):
*      Test SORCSD, SORCS2
*
*  QSV (Quotient Singular Value Decomposition):
*      Tests SGGQSV, SGGQRJ
*
*  PSV (Product Singular Value Decomposition):
*      Tests SGGPSV, SGGBRD
*
*  GSV (Generalized Singular Value Decomposition):
*      Tests SGGSVD, SGGSVP, STGSJA, SLAGS2, SLAPLL, and SLAPMT
*      For detail, see schkee.f, this test is for comparison with 
*      QSV only!
*
*  Each test path has a different set of inputs.
*  The first line of input should contain one of the
*  3-character path names in columns 1-3.  The number of remaining lines
*  depend on what is found on the first line.
*
*  The matrices used in testing are all randomly generated.
*
*-----------------------------------------------------------------------
*
*  CSD, QSV, and PSV all have the similar input file:
*
*  line 1:  'CSD' or 'QSV' or 'PSV' in columns 1 to 3.
*
*  line 2:  NN, INTEGER
*           Number of values of M, P, and N.
*
*  line 3:  MVAL, INTEGER array, dimension(NN)
*           Values of M (row dimension).
*
*  line 4:  PVAL, INTEGER array, dimension(NN)
*           Values of P (row dimension).
*
*  line 5:  NVAL, INTEGER array, dimension(NN)
*           Values of N (column dimension).
*
*  line 6:  THRESH, REAL
*           Threshold value for the test ratios.  Information will be
*           printed about each test for which the test ratio is greater
*           than or equal to the threshold.
*
*  line 7:  TSTERR, LOGICAL
*           Flag indicating whether or not to test the error exits for
*           the LAPACK routines and driver routines.
*
*  line 8:  NEWSD, INTEGER
*           A code indicating how to set the random number seed.
*           = 0:  Set the seed to a default value before each run
*           = 1:  Initialize the seed to a default value only before the
*                 first run
*           = 2:  Like 1, but use the seed values on the next line
*
*  If line 8 was 2:
*
*  line 9:  INTEGER array, dimension (4)
*           Four integer values for the random number seed.
*
*  lines 9-EOF:  Lines specifying matrix types, as for NEP.
*           The 3-character path name is 'QSV' for the generalized
*           SVD routines.      
*
*-----------------------------------------------------------------------
*
*  NMAX is currently set to 132 and must be at least 1 
*  and LWORK = NMAX*(5*NMAX+5)+1 in the parameter
*  statements below.  
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX
      PARAMETER          ( NMAX = 1320 )
      INTEGER            NCMAX
      PARAMETER          ( NCMAX = 20 )
      INTEGER            NEED
      PARAMETER          ( NEED = 14 )
      INTEGER            LWORK
      PARAMETER          ( LWORK = NMAX*( 8*NMAX+5 )+1 )
      INTEGER            LIWORK
      PARAMETER          ( LIWORK = NMAX*( 5*NMAX+20 ) )
      INTEGER            MAXIN
      PARAMETER          ( MAXIN = 20 )
      INTEGER            MAXT
      PARAMETER          ( MAXT = 30 )
      INTEGER            NIN, NOUT
      PARAMETER          ( NIN = 5, NOUT = 6 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FATAL, CSD, QSV, PSV, GSV, TSTERR
      CHARACTER          C1
      CHARACTER*3        C3, PATH
      CHARACTER*10       INTSTR
      CHARACTER*80       LINE
      INTEGER            I, I1, IC, INFO, ITMP, K, LENP, MAXTYP, NEWSD,
     $                   NK, NN, NPARMS, NRHS, NTYPES
      REAL               EPS, S1, S2, THRESH, THRSHN
*     ..
*     .. Local Arrays ..
      INTEGER            IOLDSD( 4 ), ISEED( 4 ), IWORK( LIWORK ),
     $                   MVAL( MAXIN ), NVAL( MAXIN ), PVAL( MAXIN )
      REAL               A( NMAX*NMAX, NEED ), B( NMAX*NMAX, 5 ),
     $                   C( NCMAX*NCMAX, NCMAX*NCMAX ), D( NMAX, 12 ),
     $                   RESULT( 500 ), TAUA( NMAX ), TAUB( NMAX ),
     $                   WORK( LWORK ), AB( 2*NMAX, NMAX ) 
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      REAL               SECOND, SLAMCH
      EXTERNAL           LSAMEN, SECOND, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCKCSD, SCKQSV, SCKPSV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN, MIN
*     ..
*     .. Scalars in Common
      INTEGER            NUNIT
*     ..
*     .. Data statements ..
      DATA               INTSTR / '0123456789' /
      DATA               IOLDSD / 0, 0, 0, 1 /
*     ..
*     .. Executable Statements ..
*
      S1 = SECOND( )
      FATAL = .FALSE.
      NUNIT = NOUT
*
*     Return to here to read multiple sets of data
*
   10 CONTINUE
*
*     Read the first line and set the 3-character test path
*
      READ( NIN, FMT = '(A80)', END = 380 )LINE
      PATH = LINE( 1: 3 )
      CSD = LSAMEN( 3, PATH, 'CSD' )
      QSV = LSAMEN( 3, PATH, 'QSV' )
      PSV = LSAMEN( 3, PATH, 'PSV' )
      GSV = LSAMEN( 3, PATH, 'GSV' )
*
*     Report values of parameters.
*
      IF( PATH.EQ.'   ' ) THEN
         GO TO 10
      ELSE IF( CSD ) THEN
         WRITE( NOUT, FMT = 9961 )
      ELSE IF( QSV ) THEN
         WRITE( NOUT, FMT = 9960 )
      ELSE IF( PSV ) THEN
         WRITE( NOUT, FMT = 9959 )
      ELSE IF( GSV ) THEN
         WRITE( NOUT, FMT = 9969 )
      END IF
*
*     Read the number of values of M, P, and N.
*
      READ( NIN, FMT = * )NN
      IF( NN.LT.0 ) THEN
         WRITE( NOUT, FMT = 9989 )'   NN ', NN, 1
         NN = 0
         FATAL = .TRUE.
      ELSE IF( NN.GT.MAXIN ) THEN
         WRITE( NOUT, FMT = 9988 )'   NN ', NN, MAXIN
         NN = 0
         FATAL = .TRUE.
      END IF

*     test input: 
      PRINT *, "NN = ", NN

*
*     Read the values of M
*
      READ( NIN, FMT = * )( MVAL( I ), I = 1, NN )
      DO 20 I = 1, NN
         IF( MVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9989 )' M  ', MVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( MVAL( I ).GT.NMAX ) THEN
            WRITE( NOUT, FMT = 9988 )' M  ', MVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
   20 CONTINUE
      WRITE( NOUT, FMT = 9983 )'M:    ', ( MVAL( I ), I = 1, NN )
*
*     Read the values of P
*
      READ( NIN, FMT = * )( PVAL( I ), I = 1, NN )
      DO 30 I = 1, NN
         IF( PVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9989 )' P  ', PVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( PVAL( I ).GT.NMAX ) THEN
            WRITE( NOUT, FMT = 9988 )' P  ', PVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
   30 CONTINUE
      WRITE( NOUT, FMT = 9983 )'P:    ', ( PVAL( I ), I = 1, NN )
*
*     Read the values of N
*
      READ( NIN, FMT = * )( NVAL( I ), I = 1, NN )
      DO 40 I = 1, NN
         IF( NVAL( I ).LT.0 ) THEN
            WRITE( NOUT, FMT = 9989 )' N  ', NVAL( I ), 0
            FATAL = .TRUE.
         ELSE IF( NVAL( I ).GT.NMAX ) THEN
            WRITE( NOUT, FMT = 9988 )' N  ', NVAL( I ), NMAX
            FATAL = .TRUE.
         END IF
   40 CONTINUE
      WRITE( NOUT, FMT = 9983 )'N:    ', ( NVAL( I ), I = 1, NN )
*
*     Calculate and print the machine dependent constants.
*
      WRITE( NOUT, FMT = * )
      EPS = SLAMCH( 'Underflow threshold' )
      WRITE( NOUT, FMT = 9981 )'underflow', EPS
      EPS = SLAMCH( 'Overflow threshold' )
      WRITE( NOUT, FMT = 9981 )'overflow ', EPS
      EPS = SLAMCH( 'Epsilon' )
      WRITE( NOUT, FMT = 9981 )'precision', EPS
*
*     Read the threshold value for the test ratios.
*
      READ( NIN, FMT = * )THRESH
      WRITE( NOUT, FMT = 9982 )THRESH
*
*     Read the flag that indicates whether to test the error exits.
*
      READ( NIN, FMT = * )TSTERR
*     test 
      PRINT *, "TSTERR: ", TSTERR

*
*     Read the code describing how to set the random number seed.
*
      READ( NIN, FMT = * )NEWSD
*     test
      PRINT *, "How to Set SEED: ", NEWSD
*
*     If NEWSD = 2, read another line with 4 integers for the seed.
*
      IF( NEWSD.EQ.2 ) THEN
        READ( NIN, FMT = * )( IOLDSD( I ), I = 1, 4 )
*
*		test
        print *, "New SEED:  "
        WRITE( NOUT, FMT = * )( IOLDSD( I ), I = 1, 4 )
      END IF
*
      DO 180 I = 1, 4
         ISEED( I ) = IOLDSD( I )
  180 CONTINUE
*
      IF( FATAL ) THEN
         WRITE( NOUT, FMT = 9999 )
         STOP
      END IF
*
*     Read the input lines indicating the test path and its parameters.
*     The first three characters indicate the test path, and the number
*     of test matrix types must be the first nonblank item in columns
*     4-80.
*
  190 CONTINUE
*
  200 CONTINUE
      READ( NIN, FMT = '(A80)', END = 380 )LINE
      C3 = LINE( 1: 3 )
*
      LENP = LEN( LINE )
      I = 3
      ITMP = 0
      I1 = 0
  210 CONTINUE
      I = I + 1
      IF( I.GT.LENP ) THEN
         IF( I1.GT.0 ) THEN
            GO TO 240
         ELSE
            NTYPES = MAXT
            GO TO 240
         END IF
      END IF
      IF( LINE( I: I ).NE.' ' .AND. LINE( I: I ).NE.',' ) THEN
         I1 = I
         C1 = LINE( I1: I1 )
*
*     Check that a valid integer was read
*
         DO 220 K = 1, 10
*     test
            IF( C1.EQ.INTSTR( K: K ) ) THEN
               IC = K - 1
               GO TO 230
            END IF
  220    CONTINUE
         WRITE( NOUT, FMT = 9991 )I, LINE
         GO TO 200
  230    CONTINUE
         ITMP = 10*ITMP + IC
         GO TO 210
      ELSE IF( I1.GT.0 ) THEN
         GO TO 240
      ELSE
         GO TO 210
      END IF
  240 CONTINUE
   
*     test
      NTYPES = ITMP
      PRINT *, "NTYPES = ", NTYPES

*
*     Skip the tests if NTYPES is <= 0.

      IF( NTYPES.LE.0 ) THEN
         WRITE( NOUT, FMT = 9990 )C3
         GO TO 200
      END IF
*
*     Reset the random number seed.
*
      IF( NEWSD.EQ.0 ) THEN
         DO 250 K = 1, 4
            ISEED( K ) = IOLDSD( K )
  250    CONTINUE
      END IF
*
*    =================================================
*
      IF( LSAMEN( 3, C3, 'CSD' ) ) THEN
*
*        ----------------------------------------------
*        CSD:  Cosine-Sine Decomposition
*        ----------------------------------------------
*
*         IF( TSTERR )
*     $      CALL SERRGG( 'CSD', NOUT )
          print *, "call SCKCSD ..."
         CALL SCKCSD( NN, MVAL, PVAL, NVAL, ISEED, THRESH, NMAX,
     $                AB, A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ),
     $                A( 1, 3 ), B( 1, 3 ), A( 1, 4 ), TAUA, TAUB,
     $                WORK, D( 1, 1 ), NIN, NOUT, INFO )
         IF( INFO.NE.0 )
     $      WRITE( NOUT, FMT = 9980 )'SCKCSD', INFO
*
*
      ELSE IF( LSAMEN( 3, C3, 'QSV' ) ) THEN
*
*        ----------------------------------------------
*        QSV:  Quotient Singular Value Decomposition
*        ----------------------------------------------
*
*         IF( TSTERR )
*     $      CALL SERRGG( 'QSV', NOUT )
          print *, "call SCKQSV"
         CALL SCKQSV( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX,
     $                A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ),
     $                A( 1, 3 ), B( 1, 3 ), A( 1, 4 ), TAUA, TAUB,
     $                B( 1, 4 ), IWORK, WORK, D( 1, 1 ), NIN, NOUT,
     $                INFO )
         IF( INFO.NE.0 )
     $      WRITE( NOUT, FMT = 9980 )'SCKQSV', INFO
*
*
      ELSE IF( LSAMEN( 3, C3, 'PSV' ) ) THEN
*
*        ----------------------------------------------
*        PSV:  Product Singular Value Decomposition
*        ----------------------------------------------
*
*         IF( TSTERR )
*     $      CALL SERRGG( 'PSV', NOUT )
          print *, "call SCKPSV, ........... "
         CALL SCKPSV( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX,
     $                A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ),
     $                A( 1, 3 ), B( 1, 3 ), B( 1, 4 ), TAUA,
     $                WORK, A( 1, 4 ), NIN, NOUT, INFO )
         IF( INFO.NE.0 )
     $      WRITE( NOUT, FMT = 9980 )'SCKGSV', INFO
*
      ELSE IF( LSAMEN( 3, C3, 'GSV' ) ) THEN
*
*        ----------------------------------------------
*        GSV:  Generalized Singular Value Decomposition
*        ----------------------------------------------
*
*         IF( TSTERR )
*     $      CALL SERRGG( 'GSV', NOUT )
         CALL SCKGSV( NN, MVAL, PVAL, NVAL, NTYPES, ISEED, THRESH, NMAX,
     $                A( 1, 1 ), A( 1, 2 ), B( 1, 1 ), B( 1, 2 ),
     $                A( 1, 3 ), B( 1, 3 ), A( 1, 4 ), TAUA, TAUB,
     $                B( 1, 4 ), IWORK, WORK, D( 1, 1 ), NIN, NOUT,
     $                INFO )
         IF( INFO.NE.0 )
     $      WRITE( NOUT, FMT = 9980 )'SCKGSV', INFO
*
      ELSE
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9992 )C3
      END IF
      GO TO 190
  380 CONTINUE
      WRITE( NOUT, FMT = 9994 )
      S2 = SECOND( )
      WRITE( NOUT, FMT = 9993 )S2 - S1
*
 9999 FORMAT( / ' Execution not attempted due to input errors' )
 9994 FORMAT( / / ' End of tests' )
 9993 FORMAT( ' Total time used = ', F12.2, ' seconds', / )
 9992 FORMAT( 1X, A3, ':  Unrecognized path name' )
 9991 FORMAT( / / ' *** Invalid integer value in column ', I2,
     $      ' of input', ' line:', / A79 )
 9990 FORMAT( / / 1X, A3, ' routines were not tested' )
 9989 FORMAT( ' Invalid input value: ', A6, '=', I6, '; must be >=',
     $      I6 )
 9988 FORMAT( ' Invalid input value: ', A6, '=', I6, '; must be <=',
     $      I6 )
 9983 FORMAT( 4X, A6, 10I6, / 10X, 10I6 )
 9982 FORMAT( / ' Routines pass computational tests if test ratio is ',
     $      'less than', F8.2, / )
 9981 FORMAT( ' Relative machine ', A, ' is taken to be', E16.6 )
 9980 FORMAT( ' *** Error code from ', A6, ' = ', I4 )
 9961 FORMAT( / ' Tests of the Cosine-Sine Decomposition routine,',
     $      'SORCSD. ' )
 9960 FORMAT( / ' Tests of the Quotient Singular Value Decomposition ',
     $      'routine, SGGQSV.' )
 9959 FORMAT( / ' Tests of the Product Singular Value Decomposition  ',
     $      'routine, SGGPSV. ' )
 9969 FORMAT( / ' Tests of the Generalized Singular Value',
     $      ' Decomposition routines' )
*
*     End of SCHKEE2
*
      END
