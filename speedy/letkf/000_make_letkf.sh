#!/bin/sh
set -ex
PGM=letkf020.m01
F90=$HOME/local/bin/mpif90
#JSS2#F90=f90sjx
F77=$HOME/local/bin/mpif77
OMP=
#F90OPT='-byteswapio -tp p7-64 -fast -O3'
#INLINE="-Minline"
#F90OPT='-fmax-stack-var-size=32 -fdec -fconvert=big-endian -Ofast -march=native -ffree-line-length-512 -fallow-argument-mismatch'
F90OPT='-convert big_endian -assume byterecl -O3  -march=core-avx2' # -xhost -pg'
#F90OPT='-convert big_endian -assume byterecl -O3 -check bounds'
#JSS2#F90OPT="-Umpi -O3 -Kparallel -Kdynamic_iteration -Cpp -Kprefetch_cache_level=all,prefetch_iteration_L2=50 -Ksimd -Knomfunc -Qi -Qt -Kfed"
INLINE=""
BLAS=0 #0: no blas 1: using blas

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

#$F90 $OMP $F90OPT -c spe_subfft_fftpack3.f90
$F90 $OMP $F90OPT $INLINE -c SFMT.f90
$F90 $OMP $F90OPT $INLINE -c common.f90
$F90 $OMP $F90OPT -c common_mpi.f90
$F90 $OMP $F90OPT -c common_time.f90
$F90 $OMP $F90OPT $INLINE -c common_mtx.f90
$F90 $OMP $F90OPT $INLINE -c netlib2.f
$F90 $OMP $F90OPT -c common_letkf.f90
$F90 $OMP $F90OPT -c common_lpf.f90
$F90 $OMP $F90OPT $INLINE -c common_speedy.f90
$F90 $OMP $F90OPT -c common_obs_speedy.f90
$F90 $OMP $F90OPT -c common_mpi_speedy.f90
$F90 $OMP $F90OPT -c letkf_obs.f90

$F90 $OMP $F90OPT -c $LBLAS  interpolate.f90
#$F90 $OMP $F90OPT -c filter.f90

$F90 $OMP $F90OPT -c letkf_tools.f90
$F90 $OMP $F90OPT -c letkf.f90
$F90 $OMP $F90OPT -o ${PGM} *.o #$LBLAS 


rm -f *.mod
rm -f *.o
rm -f netlib2.f
sh ulnkcommon.sh

echo "NORMAL END"
