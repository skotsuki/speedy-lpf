* help is in the end of this script
*
function markplot(args)

  if( args = '' )
    help()
    return
  endif

***** Default value *****
mark='3'
size='0.15'
color='1'
long='none'
lati='none'

***** Arguement *****
  i = 1
  while(1)
    arg = subwrd( args, i )
    i = i + 1;
    if( arg = '' ); break; endif

    while(1)
*** option
      if( arg = '-m' )  ; mark=subwrd(args,i); i=i+1; break; endif
      if( arg = '-s' )  ; size=subwrd(args,i); i=i+1; break; endif
      if( arg = '-c' )  ; color=subwrd(args,i); i=i+1; break; endif

*** lon, lat
      if( valnum(arg)!=0 & long!='none' & lati='none' ); lati=arg; break; endif
      if( valnum(arg)!=0 & long='none' ); long=arg; break; endif

      say 'syntax error : 'arg
      return

    endwhile
  endwhile

say 'mark='mark
say 'size='size
say 'color='color

'q w2xy 'long' 'lati
x1=subwrd(result,3)
y1=subwrd(result,6)
say 'long='long' -> x1='x1' lati='lati' -> y1='y1
'set line 'color
'draw mark 'mark' 'x1' 'y1' 'size

*******
* help
*
function help()
  say ' Name:'
  say '   markplot - draw mark on (lon,lat) point'
  say '               written by Y.Kamae 09/10/16'
  say ''
  say ' Usage:'
  say '   markplot long lati'
  say '           [-m mark]'
  say '           [-s size]'
  say '           [-c color]'
  say ''
  say '     mark            : style of mark (default=3:fill circle)'
  say '     size            : size of mark (default=0.15)'
  say '     color           : color of mark (default=1:black/white)'
  say ''

return
