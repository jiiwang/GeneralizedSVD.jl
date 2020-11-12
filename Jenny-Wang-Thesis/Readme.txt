Jenny Wang
University of California, Davis
December 2, 2005

This folder contains all the work of my master study in computer science at UCDavis.

JennyWang_Thesis_Final.pdf is the final version of my master thesis in pdf format. 

The software implementation of the new LAPACK-like subroutines are located in 
software/new/ and other dependent subroutines from LAPACK or BLAS are 
located in software/blas and software/lapack, respectively.

The software was tested on csia64-3.cs.ucdavis.edu
( 1.0 Ghz Itanium2 64-bit processor with 2GB of RAM ).

+ Software/  
    - Makefile  
    + new/ 
      + csd/ ( subroutines related to the CSD )
        - sorcs1.f
        - sorcs2.f
        - sorcsd.f
      + psv/ ( subroutines related to the PSVD )
        - sggbrd.f
        - sggpsv.f
      + qsv/ ( subroutines related to the QSVD )
        - sggqrj.f
        - sggpsv.f
      + tests/ ( subroutines related to testing of the CSD, PSVD, or QSVD)
        - ptmatrix.f
        - schkee2.f
        - sckcsd.f
        - sckpsv.f
        - sckqsv.f
        - scsdts.f
        - spsvts.f
        - sqsvts.f
  +  input/ ( input files for testing the CSD, PSVD, or QSVD )
  +   blas/ ( dependent subroutines from BLAS )
  + lapack/ ( dependent subroutines from LAPACK)

      

