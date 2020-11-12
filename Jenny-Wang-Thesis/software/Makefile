##----------------------------------------------------------
## UC Davis
## 2003 Jenny Wang
## Purpose: 
##   This file generates the executable for testing 
##   new LAPACK-like routines: SORCSD, SGGQSV, SGGPSV, 
##   and LAPACK routine, SGGSVD, as a reference to SGGQSV.
##
## Nov. 28, 2005
## 
## Usage: make {imkl/netlib} < infile
##----------------------------------------------------------

help:
	@echo
	@echo "Usage: make {imkl/netlib} < infile"
	@echo
	@echo "imkl   - use the Intel's Math Kernel Library 7.2.1"
	@echo "netlib - use the source codes obtained from Netlib"
	@echo "infile - see .in files in \"input\" folder"
	@echo

imkl:	use_imkl

netlib: use_netlib
	
m = lapack/modified
u = lapack/not-in-IMKL
n = new

FC = ifort -w -D_Linux 
G2C = /usr/lib/gcc-lib/ia64-redhat-linux/3.2.3
       
# modified LAPACK routines
OBJM = $m/alahdg.o $m/slatb9.o \
       $m/sgsvts.o $m/sckgsv.o $m/sggsvd.o 

# un-modified LAPACK or BLAS routines that are not in IMKL
OBJU = $u/alareq.o $u/alasum.o $u/lsamen.o $u/sbdt01.o $u/second.o \
       $u/slabad.o $u/slacpy.o $u/slagge.o $u/slagsy.o $u/slamch.o \
       $u/slaran.o $u/slarnd.o $u/slaror.o $u/slarot.o $u/slatm1.o \
       $u/slatms.o $u/sort01.o $u/ssdot.o     

# new LAPACK-style routines
OBJN = $n/csd/sorcs1.o $n/csd/sorcs2.o $n/csd/sorcsd.o \
       $n/qsv/sggqrj.o $n/qsv/sggqsv.o                 \
       $n/psv/sggbrd.o $n/psv/sggpsv.o                 \
       $n/tests/scsdts.o $n/tests/sckcsd.o             \
       $n/tests/sqsvts.o $n/tests/sckqsv.o             \
       $n/tests/spsvts.o $n/tests/sckpsv.o             \
       $n/tests/schkee2.o                              \
       $n/tests/ptmatrix.o	  
       
# ---- compile, link and execute using IMKL ----
use_imkl: $(OBJM) $(OBJU) $(OBJN)
	$(FC) $(OBJM) $(OBJU) $(OBJN) -L. \
	-L$(IMKL)/lib/64 -lmkl_lapack -lmkl_ipf -lguide  -o $@.out
	./$@.out

#-------------------------------------------------

b = blas
# dependent BLAS routines for CSD
OBJBc = $b/for-csd/blasRoutines.o    $b/for-csd/scopy.o   \
        $b/for-csd/second.o          $b/for-csd/sgemm.o   \
        $b/for-csd/snrm2.o           $b/for-csd/srot.o    \
        $b/for-csd/sscal.o           $b/for-csd/sswap.o   \
        $b/for-csd/strmm.o           $b/for-csd/xerbla.o  
	
# dependent LAPACK routines for CSD
c = lapack/for-csd
OBJC =  $c/ql/sgeql2.o $c/ql/sgeqlf.o $c/ql/sorg2l.o      \
        $c/ql/sorgql.o $c/ql/sorm2l.o $c/ql/sormql.o      \
        $c/svd/util/ieeeck.o     $c/svd/util/ilaenv.o     \
        $c/svd/single/sbdsqr.o   $c/svd/single/sgebd2.o   \
        $c/svd/single/sgebrd.o   $c/svd/single/sgelq2.o   \
        $c/svd/single/sgelqf.o   $c/svd/single/sgeqr2.o   \
        $c/svd/single/sgeqrf.o   $c/svd/single/sgesvd.o   \
        $c/svd/single/slabrd.o   $c/svd/single/slacpy.o   \
        $c/svd/single/slamch.o   $c/svd/single/slange.o   \
        $c/svd/single/slapy2.o   $c/svd/single/slarf.o    \
        $c/svd/single/slarfb.o   $c/svd/single/slarfg.o   \
	$c/svd/single/slarft.o   $c/svd/single/slartg.o   \
        $c/svd/single/slas2.o    $c/svd/single/slascl.o   \
        $c/svd/single/slaset.o   $c/svd/single/slasq1.o   \
        $c/svd/single/slasq2.o   $c/svd/single/slasq3.o   \
        $c/svd/single/slasq4.o   $c/svd/single/slasq5.o   \
        $c/svd/single/slasq6.o   $c/svd/single/slasr.o    \
        $c/svd/single/slasrt.o   $c/svd/single/slassq.o   \
        $c/svd/single/slasv2.o   $c/svd/single/sorg2r.o   \
        $c/svd/single/sorgbr.o   $c/svd/single/sorgl2.o   \
        $c/svd/single/sorglq.o   $c/svd/single/sorgqr.o   \
        $c/svd/single/sorm2r.o   $c/svd/single/sormbr.o   \
        $c/svd/single/sorml2.o   $c/svd/single/sormlq.o   \
        $c/svd/single/sormqr.o                            \
        \
	$c/tests/sasum.o         $c/tests/sbdt01.o        \
        $c/tests/slagge.o        $c/tests/slansy.o        \
        $c/tests/slaran.o        $c/tests/slarnd.o        \
        $c/tests/slarnv.o        $c/tests/slaror.o        \
        $c/tests/slaruv.o        $c/tests/sort01.o        \
        $c/tests/ssyrk.o       

# dependent LAPACK routines for GSVD
g = lapack/for-gsv
OBJG = $g/slags2.o       $g/slapll.o       $g/stgsja.o

# dependent BLAS routines for PSVD
OBJBp = $b/for-psv/sdot.o            $b/for-psv/sgemm.o   \
	$b/for-psv/snrm2.o           $b/for-psv/srot.o    \
        $b/for-psv/sscal.o           $b/for-psv/ssdot.o   \
        $b/for-psv/sswap.o           $b/for-psv/strmm.o   \
        $b/for-psv/xerbla.o

# dependent LAPACK routines for PSVD
p = lapack/for-psv
OBJP = $p/sbdsqr.o       $p/slacpy.o       $p/slamch.o   \
       $p/slapy2.o       $p/slarf.o        $p/slarfb.o   \
       $p/slarfg.o       $p/slarft.o       $p/slargv.o   \
       $p/slartg.o       $p/slas2.o        $p/slascl.o   \
       $p/slaset.o       $p/slasq1.o       $p/slasq2.o   \
       $p/slasq3.o       $p/slasq4.o       $p/slasq5.o   \
       $p/slasq6.o       $p/slasr.o        $p/slasrt.o   \
       $p/slasv2.o       $p/sorg2r.o       $p/sorgbr.o   \
       $p/sorgl2.o       $p/sorglq.o       $p/sorgqr.o   \
       \
       $p/tests/sasum.o        $p/tests/sbdt01.o         \
       $p/tests/slagge.o       $p/tests/slange.o         \
       $p/tests/slansy.o       $p/tests/slaran.o         \
       $p/tests/slarnd.o       $p/tests/slarnv.o         \
       $p/tests/slaruv.o       $p/tests/slassq.o         \
       $p/tests/sort01.o       $p/tests/ssyrk.o      

# dependent BLAS routines for QSVD
OBJBq = $b/for-qsv/blasRoutines.o    $b/for-qsv/isamax.o  \
        $b/for-qsv/scopy.o           $b/for-qsv/sgemm.o   \
        $b/for-qsv/snrm2.o           $b/for-qsv/srot.o    \
        $b/for-qsv/sscal.o           $b/for-qsv/sswap.o   \
        $b/for-qsv/strmm.o           $b/for-qsv/xerbla.o  
              
# dependent LAPACK routines for QSVD
q = lapack/for-qsv
OBJQ = $q/sgeqpf.o       $q/sgerq2.o       $q/sgerqf.o  \
       $q/sggsvp.o       $q/slapmt.o       $q/sormr2.o
	
# dependent BLAS routines 
OBJB = 	$b/blasRoutines.o $b/xerbla.o $b/sscal.o $b/srot.o  \
        $b/scopy.o        $b/sdot.o   $b/snrm2.o $b/saxpy.o \
	
# ---- compile, link and execute ----
use_netlib: $(OBJM) $(OBJU) $(OBJN) $(OBJB)   \
       libcsd libgsv libpsv libqsv
	$(FC) $(OBJM) $(OBJU) $(OBJN) $(OBJB) \
        -L.  -lcsd_lib -lqsv_lib  -lpsv_lib -lgsv_lib \
        -L$(G2C) -lg2c -o $@.out 
	./$@.out

libcsd: $(OBJBc) $(OBJC)
	$(FC) -c $(OBJBc) $(OBJC) -L$(G2C) -lg2c 
	ar rcv libcsd_lib.a $b/for-csd/*.o $c/ql/*.o \
	$c/svd/single/*.o $c/svd/util/*.o $c/tests/*.o
	ranlib libcsd_lib.a
	
libgsv: $(OBJG)
	$(FC) -c $(OBJG) -L$(G2C) -lg2c 
	ar rcv libgsv_lib.a $g/*.o
	ranlib libgsv_lib.a
	
libpsv: $(OBJBp) $(OBJP)
	$(FC) -c $(OBJBp) $(OBJBp) $(OBJP) -L$(G2C) -lg2c 
	ar rcv libpsv_lib.a $b/for-psv/*.o $p/*.o $p/tests/*.o 
	ranlib libpsv_lib.a

libqsv: $(OBJBq) $(OBJQ)
	$(FC) -c $(OBJBq) $(OBJQ) -L$(G2C) -lg2c  
	ar rcv libqsv_lib.a $b/for-qsv/*.o $q/*.o 
	ranlib libqsv_lib.a

#-------------------------------------------------
	
clean:
	find . -name '*.o' -exec rm -f {} \; 
	rm -f *.a *.out
#-------------------------------------------------

