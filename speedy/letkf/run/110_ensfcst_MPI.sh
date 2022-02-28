#!/bin/bash
#=======================================================================
# ensfcst.sh
#   This script runs the SPEEDY model with subdirectory $NODE
#=======================================================================
set -e
### input for this shell
YMDH=$1
MEMBER=$2
NODE=$3
###
if test 5$3 -eq 5 ; then
  echo "ERROR in ensfcst.sh"
  exit
fi
### run

echo $MEMBER > MEMINFO
FORT2=1111
echo $FORT2 > fort.2
echo $YMDH | cut -c1-4 >> fort.2
echo $YMDH | cut -c5-6 >> fort.2
echo $YMDH | cut -c7-8 >> fort.2
echo $YMDH | cut -c9-10 >> fort.2

#mpirun -np $NODE ./imp.exe
$MPIRUN -np $NODE ./imp.exe
#fipp -C -Ihwm -d prof mpiexec -np $NODE ./imp.exe

exit 0
