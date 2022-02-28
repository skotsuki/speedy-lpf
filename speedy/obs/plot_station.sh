#/bin/bash


for OBSNET in "raob" "reg3" "reg4" "reg5" "reg3_void" ; do
  ln -s station_${OBSNET}.tbl station.tbl.tmp

  #==> txt make
  ifort plot_station.f90 -o plot_station.exe -assume byterecl #-convert big_endian
  ./plot_station.exe

  #=== GMT ===#
  psout="./obs_map.ps"
  pngout="./obs_map.png"
  pscoast -R0/360/-90/90 -Jq180/0.028 -Ba60g30/a30g30WSne -A7000 -W1 -Di -X0.7 -K  > $psout
  psxy ./obs_map.txt -R -Jq -Sx0.1 -W3 -G0/0/0 -O >> $psout

  convert -rotate 90 -density 720 -resize 900 $psout $pngout
  convert $pngout -background white -flatten -alpha off $pngout

  mv $pngout ./figure/obs_map_${OBSNET}.png
  rm $psout 

  unlink station.tbl.tmp
done


exit
