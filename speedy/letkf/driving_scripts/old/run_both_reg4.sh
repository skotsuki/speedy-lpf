#!/bin/sh
#PBS -m abe
#PBS -M andrew.pensoneault@riken.jp
#PBS -N REG4_Interpolation_Test
module load intel/2017.4.056
module load mpt/2.16
set -e
CDIR=$PBS_O_WORKDIR
cd $CDIR/..
if [ ! -f letkf020.m01 ]; then sh make_letkf.sh || { echo 'Failed Ending' ; exit 1; }; fi;
cd run
OBS=reg4
sh 101_run_cycle_sk_sigma.sh $OBS 
