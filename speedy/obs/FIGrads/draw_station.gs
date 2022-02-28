'reinit'

inet=1 ; while(inet<=2)
  if(inet=1) ; cnet='reg2' ; endif
  if(inet=2) ; cnet='raob' ; endif
  if(inet=3) ; cnet='reg4' ; endif
  if(inet=4) ; cnet='adpu' ; endif

'set display color white'
'open Mask_project.ctl'
'reset'

'set strsiz 0.20 0.22'
if(inet=1) ; 'draw string 0.75 7.15 (a) Observing Network :: REG2 (#stations=1008)' ; endif
if(inet=2) ; 'draw string 0.75 7.15 (b) Observing Network :: RAOB (#stations=415)'  ; endif

if(inet=3) ; 'draw string 0.75 7.15 (b) Observing Network :: REG4 (#obs=264)'  ; endif
if(inet=4) ; 'draw string 0.75 7.15 (a) Observing Network :: ADPUPA (#obs=449)'  ; endif

'set lon   0 360' ; 'set xlint 30'
'set lat -90  90' ; 'set ylint 30'
*'set lon    0 360' ; 'set xlint 20'
*'set lat  -90  90' ; 'set ylint 20'

'set grid off'

'set parea 0.4 10.90 1.0 7.0'
'set grads off'
'set gxout grfill'

'./gsfile/presen.gs -alabel 0.15'

'set rgb 99 250 250 250'

'set ccols 0 99 99' 
'set clevs  0.5 1.5'
'd mask'
'./markplot_gridpnt_'cnet'.gs'
'set line 1 1 9'
'./markplot_station_'cnet'.gs'

'printim figure/station_'cnet'.png'


'!convert -gravity center -crop 780x500-5+10  figure/station_'cnet'.png  figure/station_'cnet'_gravity.png'

*'gxprint PNGfile/rmsd_'var'_'guan'_gravity.eps'


inet=inet+1 ; endwhile
return

