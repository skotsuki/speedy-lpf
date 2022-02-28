# speedy-lpf
Local Particle Filter for Numerical Weather Prediction

SPEEDY-LPF

# SPEEDY-LPF

## Quick setup
### 0. requirements

- intel fortran (ifort) or gfortran-10 (or higher version) 
- OpenMPI or MPICH for mpi
- bash
- cmake

### 1. Generate nature run
```
$ cd speedy/model/update
$ cp makefile-gfortran10 makefile # change compiler name "gfortran10" accordingly
# at the top of this system,
$ cd speedy/model/run
$ bash 000_run_first.sh  # Spin-up the nature model.
#you can expect "NORMAL END"
$ bash 001_run_cycle.sh  # actual nature run.
$ ls ../../DATA/nature/
#1981123106.grd    1982012518.grd    1982022006.grd    1982031718.grd    1982041206.grd    1982050718.grd
#1981123106_p.grd  1982012518_p.grd  1982022006_p.grd  1982031718_p.grd  1982041206_p.grd  1982050718_p.grd
#1981123112.grd    1982012600.grd    1982022012.grd    1982031800.grd    1982041212.grd    1982050800.grd
#...
```

### 2. Make observations.
```
$ cd speedy/obs
$ bash 002_make_obsmake-gfortran10.sh #for gfortran-10. change postfix accordingly.
$ bash 003_obsmake.sh #make observations, to be stored in DATA/obs
#you can expect "NORMAL END"
$ ls ../../DATA/raob/obs/
#1982010100.dat  1982012218.dat  1982021312.dat  1982030706.dat  1982032900.dat  1982041918.dat  1982051112.dat
#1982010106.dat  1982012300.dat  1982021318.dat  1982030712.dat  1982032906.dat  1982042000.dat  1982051118.dat
#1982010112.dat  1982012306.dat  1982021400.dat  1982030718.dat  1982032912.dat  1982042006.dat  1982051200.dat
```

### 3. Run DA
prepare for DA
```
$ cd speedy/letkf
$ cmake . #detect compiler automatically. 
$ make #to generate letkf020.m01, obsope.s01, and ensemble_mean_MPI.exe
# You can use ifort for parallel compilation. Use `make -j8`. but not for gfortran-10 for now(cmake bug). use serial compilation instead.
$ cd run
$ bash 000_init_kotsuki.sh
```

run DA
```
bash 200_da_letkf.sh #run SPEEDY-LETKF
bash 250_da_lpf.sh #run SPEEDY-LPF
```

you can see the result in `DATA/raob`

### Inside of DA cycle script.
`200_da_letkf.sh` and `250_da_lpf.sh` are script to run DA cycle.
See script for detail. (TBD)

### Analyze the result
TBD
