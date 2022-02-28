function presen(args)

if(subwrd(args,1)='help')
  help()
  return
endif

aint = 0
aint2 = 0
alabel = 0.25
cint = 0
clabel = 0.20
cthick = 3
mthick = 5
grads = 0

i = 1
wrd = subwrd(args,i)
while (substr(wrd,1,1) = '-')
   if (wrd = '-aint')
      i = i + 1
      aint = subwrd(args,i)
   endif
   if (wrd = '-aint2')
      i = i + 1
      aint2 = subwrd(args,i)
   endif
   if (wrd = '-alabel')
      i = i + 1
      alabel = subwrd(args,i)
   endif
   if (wrd = '-clabel')
      i = i + 1
      clabel = subwrd(args,i)
   endif
   if (wrd = '-cthick')
      i = i + 1
      cthick = subwrd(args,i)
   endif
   if (wrd = '-mthick')
      i = i + 1
      mthick = subwrd(args,i)
   endif
   if (wrd = '-grads')
      grads = 1
   endif
   i = i + 1
   wrd = subwrd(args,i)
endwhile

if(aint!=0)
'set xlint 'aint
'set ylint 'aint2
endif
if(cint!=0)
'set cint 'cint
endif

#'set mpdset hires'
'set map 1 1 'mthick
'set xlopts 1 4 'alabel
'set ylopts 1 4 'alabel
'set clopts 1 2 'clabel
'set cthick 'cthick

if(grads = 0)
'set grads off'
endif

return


** HELP **

function help()

say ''
say 'presen.gs Ver. 0.3'
say ''
say ' presen.gs [options]'
say ' ex) presen.gs'
say ''
say ' -aint (num) : interval of axis x label. (Default:No-effect)'
say ' -aint2 (num) : interval of axis y label. (Default:No-effect)'
say ' -alabel (num) : size of axis label. (Default:0.2)'
say ' -cint (num) : interval of contour label. (Default:No-effect)'
say ' -clabel (num) : size of contour label. (Default:0.15)'
say ' -cthick (num) : thickness of contour. (Default:3)'
say ' -mthick (num) : thickness of map. (Default:3)'
say ' -grads : draw grads information. (Default:off)'
say ''

return