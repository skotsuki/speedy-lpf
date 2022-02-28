#!/bin/sh
set -e
CWD=`pwd`
#####################################################################
SMPTBL="2 3 4 5 6 7 8"
#SMPTBL="7"
#####################################################################


F90OPT='-convert big_endian -assume byterecl -O3'
if [ -e prg_bilin.exe ] ; then rm prg_bilin.exe ; fi
ifort $F90OPT prg_bilin.f90 -o prg_bilin.exe

########################################### loop str
for SMP in ${SMPTBL} ; do

  if [ -e bilin.cnf ] ; then rm bilin.cnf ; fi
  cat >> bilin.cnf << EOF
  &bilin_param
    nsmp        =  ${SMP},
    dstfile     = '$CWD/inpdata/dist_all2all.grd',
  /
EOF

  ./prg_bilin.exe

  CSMP=$( printf %02i ${SMP} )
  mv markplot_bilin.gs         ./FIGrads/markplot_bilin$CSMP.gs
  mv markplot_sta3pnl_bilin.gs ./FIGrads/markplot_sta3pnl_bilin$CSMP.gs
  mv distchck.bin      ./FIGrads/distchck_bilin$CSMP.bin
  mv interp.bin        ./output/interp_bilin$CSMP.bin
done
########################################### loop end

echo "001_make-bilin.sh"
