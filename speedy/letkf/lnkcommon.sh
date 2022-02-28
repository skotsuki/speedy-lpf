#!/bin/sh
# for making link of common sources
set -e
COMMONDIR=../../common
ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/common.f90 ./
ln -fs $COMMONDIR/common_mpi.f90 ./
ln -fs $COMMONDIR/common_mtx.f90 ./
ln -fs $COMMONDIR/common_letkf.f90 ./
ln -fs $COMMONDIR/common_lpf.f90   ./
ln -fs $COMMONDIR/netlib.f ./
ln -fs $COMMONDIR/netlibblas.f ./
ln -fs $COMMONDIR/common_time.f90 ./

ln -fs ../common/common_speedy.f90 ./
ln -fs ../common/common_mpi_speedy.f90 ./
ln -fs ../common/common_obs_speedy.f90 ./

## weight smoother :: following P. Andrew 
#PARMDIR='../model/source'
#ln -fs $PARMDIR/spe_subfft_fftpack.f  ./
#ln -fs $PARMDIR/spe_spectral.f        ./
#ln -fs $PARMDIR/spe_subfft_fftpack.f  ./

#PARMDIR='../model/tmp'
#ln -fs $PARMDIR/atparam.h ./
#ln -fs $PARMDIR/atparam1.h ./
#ln -fs $PARMDIR/com_spectral.h ./


