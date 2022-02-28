
set -e
TAG=`date '+%Y%m%d'`
VER="v1.00"
BASE=`pwd`
DIR="SYS${TAG}_SPEEDY_LPF${VER}"
CP="cp -arv"

mkdir $DIR ; cd $DIR
${CP} ../archive.sh         ./
${CP} ../common/            ./

mkdir speedy
  ${CP} ../speedy/common       ./speedy/
  ${CP} ../speedy/FIGquick     ./speedy/
  ${CP} ../speedy/letkf        ./speedy/
  ${CP} ../speedy/model        ./speedy/
  ${CP} ../speedy/obs          ./speedy/
  ${CP} ../speedy/sk_pps_bla   ./speedy/ 
  ${CP} ../speedy/src_LTERPLUT ./speedy/



echo "============================================================================"
echo "                         Archived   Successfully                            "
echo "============================================================================"

