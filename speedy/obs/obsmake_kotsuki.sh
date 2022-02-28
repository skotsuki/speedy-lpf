#!/bin/sh
CDIR=`pwd`
cd ..
SPEEDY=`pwd`

#EXP=reg3
#EXP=raob
for EXP in "raob" "reg2" "reg3" "reg3_void"  "reg4" "reg5" ; do

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
	EY=1983
	EM=01
	ED=01
	EH=00
	# start
	cd $CDIR
	mkdir -p $OBSDIR
	ln -s $SPEEDY/common/orography_t30.dat fort.21
	# main loop
	while test $IY$IM$ID$IH -le $EY$EM$ED$EH
	do
	ln -s $TRUEDIR/$IY$IM$ID$IH.grd true.grd

	ln -s station_${EXP}.tbl station.tbl # kotsuki
	./$PGM
	unlink station.tbl                 # kotsuki
	
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
	rm fort.21
done
echo "NORMAL END"
