#!/bin/sh
#=======================================================================
# init.sh
#   This script prepares for new LETKF cycle-run experiment
#=======================================================================
set -e
#-----------------------------------------------------------------------
# Modify below according to your environment
#-----------------------------------------------------------------------
MEMBER=128
M=$( printf %06i ${MEMBER} )
#OBS=raob
#EXP=M${M}L500IADP

### directory settings
cd ../..
SPEEDY=`pwd`
NATURE=$SPEEDY/DATA/nature     # nature run
#OUTPUT=$SPEEDY/DATA/$OBS/$EXP  # directory for new experiment
OUTPUT=$SPEEDY/DATA/init_ensemble
### initial date setting
IYYYY=1982
IMM=01
IDD=01
IHH=00
#-----------------------------------------------------------------------
# Usually do not modify below
#-----------------------------------------------------------------------
source $SPEEDY/../common/timeinc.sh
### clean
rm -rf $OUTPUT
### mkdir
mkdir -p $OUTPUT/log
mkdir -p $OUTPUT/infl_mul
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/infl_mul
# SK for LPF#
mkdir -p $OUTPUT/peff_lpf  ;  cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/peff_lpf/
mkdir -p $OUTPUT/rsmp_lpf  ;  cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/rsmp_lpf/
mkdir -p $OUTPUT/asis_lpf  ;  cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/asis_lpf/
mkdir -p $OUTPUT/gues/kldv ;  cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/gues/kldv
mkdir -p $OUTPUT/anal/kldv ;  cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/anal/kldv

MEM=1
while test $MEM -le $MEMBER
do
#if test $MEM -lt 100
#then
#MEM=0$MEM
#fi
#if test $MEM -lt 10
#then
#MEM=0$MEM
#fi
MEM=$( printf %06i ${MEM} )
mkdir -p $OUTPUT/anal/$MEM
mkdir -p $OUTPUT/anal_f/$MEM
mkdir -p $OUTPUT/gues/$MEM
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/anal/$MEM
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/anal_f/$MEM
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/gues/$MEM
  echo "end MKDIR MEMBER::" $MEM
MEM=`expr $MEM + 1`
done

for MEM in mean sprd
do
mkdir -p $OUTPUT/anal/$MEM
mkdir -p $OUTPUT/anal_f/$MEM
mkdir -p $OUTPUT/gues/$MEM
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/anal/$MEM
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/anal_f/$MEM
cp $SPEEDY/common/yyyymmddhh*.ctl $OUTPUT/gues/$MEM
done
### copy initial conditions
TY=1982
TM=02
TD=01
TH=00
MEM=1
while test $MEM -le $MEMBER
do
#if test $MEM -lt 100
#then
#MEM=0$MEM
#fi
#if test $MEM -lt 10
#then
#MEM=0$MEM
#fi
MEM=$( printf %06i ${MEM} )
cp $NATURE/$TY$TM$TD$TH.grd $OUTPUT/gues/$MEM/$IYYYY$IMM$IDD$IHH.grd
UY=`timeinc6hr $TY $TM $TD $TH | cut -c1-4`
UM=`timeinc6hr $TY $TM $TD $TH | cut -c5-6`
UD=`timeinc6hr $TY $TM $TD $TH | cut -c7-8`
UH=`timeinc6hr $TY $TM $TD $TH | cut -c9-10`
TY=`timeinc6hr $UY $UM $UD $UH | cut -c1-4`
TM=`timeinc6hr $UY $UM $UD $UH | cut -c5-6`
TD=`timeinc6hr $UY $UM $UD $UH | cut -c7-8`
TH=`timeinc6hr $UY $UM $UD $UH | cut -c9-10`
  echo "END CP MEM :: " $MEM "by" $NATURE/$TY$TM$TD$TH.grd
MEM=`expr $MEM + 1`
done

echo "NORMAL END"
