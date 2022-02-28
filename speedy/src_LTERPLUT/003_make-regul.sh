#!/bin/sh
set -e
CWD=`pwd`
#####################################################################
#SMPTBL="2 3 4 5 6 7 8 9"
#SMPTBL="2 3 4 5 6 7"
#SMPTBL="2 3 4 5 6 7 8"
SMPTBL="2 3 4 5 6 7 8"
#####################################################################


F90OPT='-convert big_endian -assume byterecl -O3'
if [ -e prg_regul.exe ] ; then rm prg_regul.exe ; fi
ifort $F90OPT prg_regul.f90 -o prg_regul.exe

########################################### loop str
for SMP in ${SMPTBL} ; do

  if [ -e regul.cnf ] ; then rm regul.cnf ; fi
  cat >> regul.cnf << EOF
  &regul_param
    nsmp        =  ${SMP},
    dstfile     = '$CWD/inpdata/dist_all2all.grd',
    grdfile     = '$CWD/inpdata/grid_all2all.grd',
    raobtbl     = '$CWD/../obs/station_raob.tbl'
  /
EOF

  ./prg_regul.exe

  CSMP=$( printf %02i ${SMP} )
  mv markplot_regul.gs         ./FIGrads/markplot_regul$CSMP.gs
  mv markplot_sta3pnl_regul.gs ./FIGrads/markplot_sta3pnl_regul$CSMP.gs
  mv distchck.bin      ./FIGrads/distchck_regul$CSMP.bin
  mv interp.bin        ./output/interp_regul$CSMP.bin
done
########################################### loop end

echo "003_make-regul.sh"
