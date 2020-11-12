       SUBROUTINE PTMATRIX( M, N, A, LDA )
         INTEGER M, N, I, J
         REAL    A( LDA, * )
         print *, " [ "
         DO 10 I = 1, M
              WRITE( *, 900 ) ( A( I, J ), J = 1, N )
  10     CONTINUE
         print *, " ]; "
  900    FORMAT(8F)
       RETURN
       END
