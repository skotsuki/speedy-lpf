cmake_minimum_required(VERSION 3.0)
# cmake 自体のバージョンを明示します。

enable_language(Fortran)

project(speedy-letkf)

SET(CMAKE_F90_COMPILER "mpif90")
# SET(CMAKE_F90_FFLAGS_RELEASE "-convert big_endian -assume byterecl -O3")
SET(CMAKE_F90_FFLAGS_DEBUG "-g")

include(CMakeForceCompiler)

set(CMAKE_Fortran_COMPILER mpif90)

message("compiler: " ${CMAKE_FORTRAN_COMPILER_ID})

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  SET(CMAKE_Fortran_FLAGS "-convert big_endian -assume byterecl -O3  -march=core-avx2")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-fdec -fconvert=big-endian -O3 -ffree-line-length-512 -fallow-argument-mismatch")
else()
  message(unknow compiler)
endif()
  
set(COMMON ../../common)
set(COMMON_SPEEDY ../common)

add_executable(letkf020.m01
${COMMON}/SFMT.f90
${COMMON}/common.f90
${COMMON}/common_mpi.f90
${COMMON}/common_time.f90
${COMMON}/common_mtx.f90
${COMMON}/netlib2.f
${COMMON}/common_letkf.f90
${COMMON}/common_lpf.f90
${COMMON_SPEEDY}/common_speedy.f90
${COMMON_SPEEDY}/common_obs_speedy.f90
${COMMON_SPEEDY}/common_mpi_speedy.f90
./letkf_obs.f90
./interpolate.f90
./letkf_tools.f90
./letkf.f90
)


add_executable(obsope_MPI.s01
${COMMON}/SFMT.f90
${COMMON}/common.f90
${COMMON}/common_mpi.f90
${COMMON}/common_mtx.f90
${COMMON}/netlib2.f
${COMMON}/common_letkf.f90
${COMMON_SPEEDY}/common_speedy.f90
${COMMON_SPEEDY}/common_obs_speedy.f90
./obsope_MPI.f90
)

add_executable(prg_ensemble_mean_MPI.exe
./prg_ensemble_mean_MPI.f90
)
