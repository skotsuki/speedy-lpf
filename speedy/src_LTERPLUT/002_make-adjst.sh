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
if [ -e prg_adjst.exe ] ; then rm prg_adjst.exe ; fi
ifort $F90OPT prg_adjst.f90 -o prg_adjst.exe

########################################### loop str
for SMP in ${SMPTBL} ; do

  if [ -e adjst.cnf ] ; then rm adjst.cnf ; fi
  cat >> adjst.cnf << EOF
  &adjst_param
    nsmp        =  ${SMP},
    dstfile     = '$CWD/inpdata/dist_all2all.grd',
    grdfile     = '$CWD/inpdata/grid_all2all.grd',
    raobtbl     = '$CWD/../obs/station_raob.tbl'
  /
EOF

  ./prg_adjst.exe

  CSMP=$( printf %02i ${SMP} )
  mv markplot_adjst.gs ./FIGrads/markplot_adjst$CSMP.gs
  mv distchck.bin      ./FIGrads/distchck_adjst$CSMP.bin
  mv interp.bin        ./output/interp_adjst$CSMP.bin
done
########################################### loop end

echo "002_make-adjst.sh"
