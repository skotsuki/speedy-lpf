#!/bin/bash

#======================================
ymdh=$1   ; ymdh=${ymdh:=1982010100}
member=$2 ; member=${member:=20}
EXP=$3    ; EXP=${EXP:=M100L1000IADP_MPI}
OUTPUT=$4 ; OUTPUT=${OUTPUT:=~/miyoshi-read-only/speedy/DATA/raob/$EXP}
NODE=$5     ; NODE=${NODE:=4}
#======================================

#mpifrtpx src/prg_ensemble_mean_MPI.f90 -o bin/ensemble_mean_MPI.exe -Kfast,parallel,noeval,ocl -V -Qa,d,i,p,t,x -Koptmsg=2

### gues
#data_dir="$OUTPUT/gues"
#echo $ymdh     >  info.txt
#echo $member   >> info.txt
#echo $data_dir >> info.txt
#mpiexec -np $NODE ./ensemble_mean_MPI.exe

### anal_f
data_dir="$OUTPUT/anal_f"
echo $ymdh     >  info.txt
echo $member   >> info.txt
echo $data_dir >> info.txt
#mpirun -np $NODE ./ensemble_mean_MPI.exe
$MPIRUN -np $NODE ./ensemble_mean_MPI.exe
exit
