'reinit' 
'set display color white'
'open resample.ctl'
*******************************************************************************************
'reset'
'set display color white'

'set grads off'


*******************************************************************************************
'./gsfile/define_colors.gs'
'set rgb 69 164  32  32'
'set rgb 65 255 100 200'
'set rgb 63 255 200 255'

'set lon 0.0 40.0'
'set lat 0.0 40.0'
'set xlint 5'
'set ylint 5'
'set yflip on'
'set mproj off'

'set gxout grfill'
'./gsfile/presen.gs -alabel 0.12'
*******************************************************************************************
xlll0=0.8 ; xrrr0=0.9
xlll1=1.4 ; xrrr1=5.4
xlll2=5.9 ; xrrr2=9.9
yddd1=4.5 ; yuuu1=8.0
yddd2=0.3 ; yuuu2=3.8

*******************************************************************************************
'set strsiz 0.16 0.18'
'draw string 0.5     8.12 (a) w'
'draw string 'xlll1' 8.12 (b) Sampled T (1st sample)'
'draw string 'xlll2' 8.12 (c) Sampled T (2nd sample)'
'draw string 'xlll1' 3.92 (d) Monte-Carlo T (N=200)'
'draw string 'xlll2' 3.92 (e) Monte-Carlo T (N=10000)'

*******************************************************************************************
ic=1 ; while(ic<=5)
*  'set parea 2.0 9.5 1.0 7.9'
  if(ic=1|ic=3); xlll=xlll1 ; xrrr=xrrr1 ; endif
  if(ic=2|ic=4); xlll=xlll2 ; xrrr=xrrr2 ; endif
  if(ic=1|ic=2); yddd=yddd1 ; yuuu=yuuu1 ; endif
  if(ic=3|ic=4); yddd=yddd2 ; yuuu=yuuu2 ; endif
  'set parea 'xlll' 'xrrr' 'yddd' 'yuuu''
  if(ic=5)
    'set lon 0.0      1.0'
    'set lat 0.00001 40.0'
    'set xlab off'
    'set parea 'xlll0' 'xrrr0' 'yddd1' 'yuuu1''
  endif

*  'set clevs 0.01 0.02 0.03 0.04 0.05 0.1 0.2 0.3 1.0'
*  'set ccols 0 5 11 10 7 8 2 69 6 65 63'

  'set clevs 0.001 0.01 0.02 0.03 0.04 0.05 0.1 0.2 0.3 0.4 1.0'
  'set ccols 0 41 43 5 11 10 7 8 2 69 6 65 63'


  if(ic=1);'d mat001+0.000001';endif
  if(ic=2);'d mat002+0.000001';endif
  if(ic=3);'d kkk200';endif
  if(ic=4);'d kkk999';endif
  if(ic=5);'d weight';endif

  if(ic=1) ; './gsfile/xcbar.gs 10.20 10.35  0.3 8.0 -dir v  -line on -edge triangle -fwidth 0.11 -fheight 0.12' ; endif
ic=ic+1
endwhile

'printim PNGfile/matrix.png'
return