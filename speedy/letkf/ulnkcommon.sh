#!/bin/sh
# for deleting links of common sources
set -e

rm -f SFMT.f90
rm -f common.f90
rm -f common_mpi.f90
rm -f common_mtx.f90
rm -f common_letkf.f90
rm -f common_lpf.f90
rm -f netlib.f
rm -f netlibblas.f
rm -f common_time.f90

rm -f common_speedy.f90
rm -f common_mpi_speedy.f90
rm -f common_obs_speedy.f90



#rm -f spe_subfft_fftpack.f 
#rm -f spe_subfft_fftpack2.f 



#rm -f atparam.h 
#rm -f atparam1.h 
#rm -f com_spectral.h 
