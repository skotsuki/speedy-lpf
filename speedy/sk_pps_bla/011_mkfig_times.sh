#!/bin/bash

CWD=`pwd`
rm ./fig_011/*.eps ./fig_011/*.png
gnuplot "gscript_011.gnu"

###==> convert fig.
cd ./fig_011/

FNUMTMP=`find ./ -name "*.eps" -print`
#FNUMTMP=`find ./ -name "????.eps" -print`
for FNUM in ` echo ${FNUMTMP}` ; do
#  FNAME=`echo ${FNUM}|cut -c 1-13`
  basename=${FNUM##*/}
  filename=${basename%.*}
  FOUT=${filename}.png
  FEPS=${filename}.eps
  echo "   converting... " $FNUM " ==>  " $FOUT

#  convert -density 720 -resize 900 $FNUM $FOUT
  convert -density 120 -resize 900 $FNUM $FOUT
  convert $FOUT -background white -flatten -alpha off ${FOUT}
  convert -gravity center -crop 880x510+10+00 ${FOUT} ${FOUT}
  rm $FNUM  
done
echo "   montaging..."
montage time_reg2_rmse_lpf.png            time_reg2_peff_lpf.png           -tile 1x2 -geometry 880x510 montage_reg2_lpf.png
montage time_reg2_rmse_lpfgm_fgt1.00.png  time_reg2_peff_lpfgm_fgt1.00.png -tile 1x2 -geometry 880x510 montage_reg2_lpfgm.png
montage time_raob_rmse_lpfgm_fgt1.00.png  time_raob_peff_lpfgm_fgt1.00.png -tile 1x2 -geometry 880x510 montage_raob_lpfgm.png
exit

