'reinit'
case=1 ; while(case<=2)
  if(case=1) ; cnam='bilin' ; endif
  if(case=2) ; cnam='adjst' ; endif
  if(case=3) ; cnam='regul' ; endif
inet=2 ; while(inet<=8)
  if(inet=2) ; cnet='02' ; endif
  if(inet=3) ; cnet='03' ; endif
  if(inet=4) ; cnet='04' ; endif
  if(inet=5) ; cnet='05' ; endif
  if(inet=6) ; cnet='06' ; endif
  if(inet=7) ; cnet='07' ; endif
  if(inet=8) ; cnet='08' ; endif
  if(inet=9) ; cnet='09' ; endif

  '! rm distchck_ref.bin'
  fname='distchck_'cnam''cnet'.bin'
  '! ln -s 'fname' distchck_ref.bin '

'open distchck.ctl'
'set display color white'
'reset'



'set strsiz 0.18 0.20'
if(case=1)
  if(inet=2) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: uniform, #02 (1200/4608)'  ; endif
  if(inet=3) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: uniform, #03 (0544/4608)'  ; endif
  if(inet=4) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: uniform, #04 (0312/4608)'  ; endif
  if(inet=5) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: uniform, #05 (0220/4608)'  ; endif
  if(inet=6) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: uniform, #06 (0144/4608)'  ; endif
  if(inet=7) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: uniform, #07 (0112/4608)'  ; endif
  if(inet=8) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: uniform, #08 (0084/4608)'  ; endif
  if(inet=9) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: uniform, #09 (0077/4608)'  ; endif
endif
if(case=2|case=3)
  if(inet=2) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: regulated, #02 (1200/4608)'  ; endif
  if(inet=3) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: regulated, #03 (0544/4608)'  ; endif
  if(inet=4) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: regulated, #04 (0312/4608)'  ; endif
  if(inet=5) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: regulated, #05 (0220/4608)'  ; endif
  if(inet=6) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: regulated, #06 (0144/4608)'  ; endif
  if(inet=7) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: regulated, #07 (0112/4608)'  ; endif
  if(inet=8) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: regulated, #08 (0084/4608)'  ; endif
  if(inet=9) ; 'draw string 0.75 8.15 Reference Points for Weight Interp. :: regulated, #09 (0077/4608)'  ; endif
endif


'set lon   0 360'
'set lat -90  90'
'set grid off'

'set parea 0.4 10.90 1.0 7.0'
'set grads off'
'set gxout grfill'

'./gsfile/presen.gs -alabel 0.15'

*'set lat -85 85'

'set parea 0.4 10.90 4.5 8.0' ; 'd dinp'      ; 
  'xcbar  9.0  9.15 2.50 8.00 -dir v -line on -edge triangle -fwidth 0.11 -fheight 0.11'
'set parea 0.4 10.90 0.5 4.0' ; 'd abs(dinp-dist)'
  'xcbar 10.0 10.15 0.50 6.00 -dir v -line on -edge triangle -fwidth 0.11 -fheight 0.11'

'set rgb 99 250 250 250'

*'./markplot_gridpnt.gs'
*'set line 1 1 9'
*'./markplot_bilin'cnet'.gs'

'printim figure/dinp_'cnam'_'cnet'.png'


*'!convert -gravity center -crop 780x500-5+10  figure/bilinear_'cnet'.png  figure/station_'cnet'_gravity.png'

*'gxprint PNGfile/rmsd_'var'_'guan'_gravity.eps'


inet=inet+1 ; endwhile
case=case+1 ; endwhile
return

