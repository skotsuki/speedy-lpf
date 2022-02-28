F90=$HOME/local/bin/mpif90
#JSS2#F90=f90sjx
F77=$HOME/local/bin/mpif77
F90OPT='-convert big_endian -assume byterecl -O3' # -traceback -check all' # -pg'

$F90 prg_ensemble_mean_MPI.f90 -o ensemble_mean_MPI.exe $F90OPT
#JSS2#F90OPT="-Umpi -O3"
#JSS2#f90sjx $F90OPT prg_ensemble_mean_MPI.f90 -o ensemble_mean_MPI.exe

echo "NORMAL END"
