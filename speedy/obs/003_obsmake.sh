#!/bin/sh
CDIR=`pwd`
cd ..
SPEEDY=`pwd`
for EXP in "reg2" "reg4" "raob"; do
#for EXP in "reg4" ; do
##
cd $CDIR
ln -fs station_${EXP}.tbl station.tbl
cd ../
##

TRUEDIR=$SPEEDY/DATA/nature
OBSDIR=$SPEEDY/DATA/$EXP/obs
PGM=obsmake.s01
source $SPEEDY/../common/timeinc.sh
# Initial date
IY=1982
IM=01
ID=01
IH=00
# Final date
EY=1982
EM=06
ED=01
EH=00
# start
cd $CDIR
mkdir -p $OBSDIR
#ln -fs $SPEEDY/common/orography_t30.dat fort.21
ln -fs $SPEEDY/common/orography_t30.dat orography_t30.dat
# main loop
while test $IY$IM$ID$IH -le $EY$EM$ED$EH
do
ln -sf $TRUEDIR/$IY$IM$ID$IH.grd true.grd

./$PGM

mv obs.dat $OBSDIR/$IY$IM$ID$IH.dat
rm true.grd

TY=`timeinc6hr $IY $IM $ID $IH | cut -c1-4`
TM=`timeinc6hr $IY $IM $ID $IH | cut -c5-6`
TD=`timeinc6hr $IY $IM $ID $IH | cut -c7-8`
TH=`timeinc6hr $IY $IM $ID $IH | cut -c9-10`
IY=$TY
IM=$TM
ID=$TD
IH=$TH
done
#rm fort.21
rm station.tbl
echo "NORMAL END"
done
