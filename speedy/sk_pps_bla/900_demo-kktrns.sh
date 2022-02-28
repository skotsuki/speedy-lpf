#!/bin/sh
set -e
CWD=`pwd`

# ========>>> preparation           <<<============#
FOPT='-convert big_endian -assume byterecl -O3 -debug all'
###FOPT='-convert big_endian -assume byterecl -fopenmp'
if [ -e prg_demo-kktrns.exe ] ; then rm prg_demo-kktrns.exe ; fi
#ifort $FOPT prg_demo-kktrns.f90 -o prg_demo-kktrns.exe
ifort $FOPT prg_demo-mptrns.f90 -o prg_demo-kktrns.exe

./prg_demo-kktrns.exe


echo "end 900_demo-kktrns.sh"
exit