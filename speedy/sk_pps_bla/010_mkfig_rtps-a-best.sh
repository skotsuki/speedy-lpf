#!/bin/bash

CWD=`pwd`
rm ./fig_010/*.eps ./fig_010/*.png
gnuplot "gscript_010.gnu"

###==> convert fig.
cd ./fig_010/

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
  convert -gravity center -crop 850x630+00+0 ${FOUT} ${FOUT}
  rm $FNUM  
done 
montage rtps_reg2_2LAPF__all_Hxxxx.png   rtps_reg2_2LAPF__lim_Hxxxx.png  -tile 2x1 -geometry 850x630 montage__lapf_reg2.png
montage rtps_reg2_3LPFGM_tau0.0_lim.png  rtps_reg2_3LPFGM_tau1.0_lim.png -tile 2x1 -geometry 850x630 montage_lpfgm_reg2.png
montage rtps_raob_3LPFGM_tau1.0_all.png  rtps_raob_3LPFGM_tau1.0_lim.png -tile 2x1 -geometry 850x630 montage_lpfgm_raob.png

#montage rtps_a.png rtps_b.png -tile 2x1 -geometry 10x10 montage_rtps.png


exit
