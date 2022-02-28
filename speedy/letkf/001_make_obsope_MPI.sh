#!/bin/sh
set -ex
PGM=obsope_MPI.s01
#F90=mpifrtpx
F90=$HOME/local/bin/mpif90
#JSS2#F90=f90sjx
F77=$HOME/local/bin/mpif77
#JSS2#F90=f90sjx
OMP=
#F90OPT='-byteswapio -tp sandybridge-64 -fast -O3'
#F90OPT='-fdec -fconvert=big-endian -O3 -ffree-line-length-512 -fallow-argument-mismatch'
F90OPT='-convert big_endian -assume byterecl -O3' # -traceback -check all' # -pg'
#JSS2#F90OPT="-Umpi -O3 -Kparallel -Kdynamic_iteration -Cpp -Kprefetch_cache_level=all,prefetch_iteration_L2=50 -Ksimd -Knomfunc -Qi -Qt -Kfed"
#F90OPT='-Kfast,parallel,noeval,ocl,simd -V -Qa,d,i,p,t,x'
#INLINE="-Minline"
BLAS=1 #0: no blas 1: using blas

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

cat netlib.f > netlib2.f
if test $BLAS -eq 1
then
LBLAS="-mkl"
cat netlibblas.f >> netlib2.f #remove soon
else
cat netlibblas.f >> netlib2.f
LBLAS=""
fi

$F90 $OMP $F90OPT $INLINE -c SFMT.f90
$F90 $OMP $F90OPT $INLINE -c common.f90
$F90 $OMP $F90OPT -c common_mpi.f90
$F90 $OMP $F90OPT -c common_mtx.f90
$F90 $OMP $F90OPT $INLINE -c netlib2.f
$F90 $OMP $F90OPT -c common_letkf.f90
$F90 $OMP $F90OPT $INLINE -c common_speedy.f90
$F90 $OMP $F90OPT -c common_obs_speedy.f90
$F90 $OMP $F90OPT -c obsope_MPI.f90
$F90 $OMP $F90OPT -o ${PGM} *.o

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh

echo "NORMAL END"
