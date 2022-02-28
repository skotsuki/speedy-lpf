#!/bin/bash

CWD=`pwd`
rm ./fig_910/*.eps ./fig_910/*.png
gnuplot "gscript_910.gnu"

###==> convert fig.
cd ./fig_910/

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
#  convert -gravity center -crop 750x630+00+0 ${FOUT} ${FOUT}
  rm $FNUM  
done 


exit
