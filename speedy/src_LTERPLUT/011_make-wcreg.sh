#!/bin/sh
set -e
CWD=`pwd`
#####################################################################
#SMPTBL="2 3 4 5 6 7 8 9"
#SMPTBL="2 3 4 5 6 7"
#SMPTBL="2 3 4 5 6 7 8"
SMPTBL="2 3 4 5 6"
MEMTBL="16"
##################################### fixed
OBSTBL="raob"
PERIOD="1982020100-1982022818"
WCORDIR="/home/kotsuki/daschool_chk/letkf-master_hakushu_sk6_xlng/speedy/DATA"
#####################################################################


F90OPT='-convert big_endian -assume byterecl -O3'
if [ -e prg_wcreg.exe ] ; then rm prg_wcreg.exe ; fi
ifort $F90OPT prg_wcreg.f90 -o prg_wcreg.exe

########################################### loop str
for MEM in ${MEMTBL} ; do
for OBS in ${OBSTBL} ; do
  if [ $MEM -eq 16 ] ; then
      if   [ $OBS == "reg2" ] ; then OPTHSG=500
      elif [ $OBS == "reg4" ] ; then OPTHSG=600
      elif [ $OBS == "raob" ] ; then OPTHSG=700
      else echo "  error stop, no such memeber & obs implemented :: " $MEM $OBS
      fi
    elif [ $MEM -eq 32 ] ; then
      if   [ $OBS == "reg2" ] ; then OPTHSG=500
      elif [ $OBS == "reg4" ] ; then OPTHSG=700
      elif [ $OBS == "raob" ] ; then OPTHSG=800
      else echo "  error stop, no such memeber & obs implemented :: " $MEM $OBS
      fi
    elif [ $MEM -eq 64 ] ; then
      if   [ $OBS == "reg2" ] ; then OPTHSG=600
      elif [ $OBS == "reg4" ] ; then OPTHSG=700
      elif [ $OBS == "raob" ] ; then OPTHSG=1200
      else echo "  error stop, no such memeber & obs implemented :: " $MEM $OBS
      fi
    elif [ $MEM -eq 128 ] ; then
      if   [ $OBS == "reg2" ] ; then OPTHSG=800
      elif [ $OBS == "reg4" ] ; then OPTHSG=900
      elif [ $OBS == "raob" ] ; then OPTHSG=1300
      else echo "  error stop, no such memeber & obs implemented :: " $MEM $OBS
      fi
    else
      echo "  error stop, no such memeber implemented :: " $MEM
      exit
    fi
    HSIG=$OPTHSG

    M=$( printf %06i ${MEM} )
    H=$( printf %04i ${HSIG}   )
    WCORDATA=$WCORDIR/$OBS/M${M}L${H}IADP/postanl/wcorall2all_LEV004_SMP${PERIOD}.grd

    if [ ! -e $WCORDATA ]; then 
      echo " error, no such file :: " $WCORDATA
      exit
    fi
########################### main part start
for SMP in ${SMPTBL} ; do

  if [ -e wcreg.cnf ] ; then rm wcreg.cnf ; fi
  cat >> wcreg.cnf << EOF
  &wcreg_param
    nsmp        =  ${SMP},
    dstfile     = '$CWD/inpdata/dist_all2all.grd',
    grdfile     = '$CWD/inpdata/grid_all2all.grd',
    raobtbl     = '$CWD/../obs/station_raob.tbl'
    wcrfile     = '${WCORDATA}',
  /
EOF

  ./prg_wcreg.exe

  CSMP=$( printf %02i ${SMP} )
  mv interp.bin        ./output/interp_wcreg${CSMP}_${OBS}H${H}.bin
  #mv distchck.bin      ./FIGrads/distchck_wcreg$CSMP.bin
  #mv interp.bin        ./output/interp_wcreg$CSMP.bin
done # SMP
done # MEM
done # MEM
########################################### loop end

echo "011_make-wcreg.sh"
