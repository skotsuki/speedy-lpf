*
* Help is in the end of this script.
*
function xcbar( args )
  _version='0.07r1'

  if( args = '' )
    help()
*    return
  endif

***** Default value *****

  xmin = 'none'
  xmax = 'none'
  ymin = 'none'
  ymax = 'none'
  direction = 'none'

  cnum = -1


***** Arguement *****

* position
*  xmin = subwrd( args, 1 )
*  xmax = subwrd( args, 2 )
*  ymin = subwrd( args, 3 )
*  ymax = subwrd( args, 4 )
*  temp = subwrd( args, 5 )

*  if( valnum(temp) = 0 )
*** Usage 1 ***
*    i = 5
    i = 1
    arg = 'dummy'
    while( 1 )

      arg = subwrd( args, i )
      i = i + 1
      if( arg = '' ); break; endif

      while( 1 )

        if( arg = '-direction' | arg = '-dir' )
          direction = subwrd( args, i )
          if( direction = 'h' ); direction = 'horizontal'; endif
          if( direction = 'v' ); direction = 'vertical'; endif
          i = i + 1
          break
        endif
        if( arg = '-edge' )
          edge = subwrd( args, i )
          i = i + 1
          break
        endif
        if( arg = '-fwidth' | arg = '-fw' )
          fwidth = subwrd( args, i )
          i = i + 1
          break
        endif
        if( arg = '-fheight' | arg = '-fh' )
          fheight = subwrd( args, i )
          i = i + 1
          break
        endif
        if( arg = '-fthickness' | arg = '-ft' )
          fthickness = subwrd( args, i )
          i = i + 1
          break
        endif
        if( arg = '-fstep' | arg = '-fs' )
          fstep = subwrd( args, i )
          i = i + 1
          break
        endif
        if( arg = '-foffset' | arg = '-fo' )
          foffset = subwrd( args, i )
          i = i + 1
          break
        endif
        if( arg = '-line' )
          line = subwrd( args, i )
          i = i + 1
          break
        endif

        if( arg = '-levcol' )
*         color(1) level(1) color(2) level(2) ... level(cnum-1) color(cnum)
          cnum = 1
          while( 1 )
            col.cnum = subwrd( args, i )
            if( valnum(col.cnum) = 0 )
              say 'Error in xcbar.gs: Syntax error in -levcol'
              return
            endif
            i = i + 1
            hi.cnum = subwrd( args, i )
            if( valnum(hi.cnum) = 0 ); break; endif
            cnum = cnum + 1
            i = i + 1
          endwhile
          break
        endif

        if( valnum(arg) != 0 & xmin = 'none' )
          xmin = arg
          break
        endif

        if( valnum(arg) != 0 & xmax = 'none' )
          xmax = arg
          break
        endif
        
        if( valnum(arg) != 0 & ymin = 'none' )
          ymin = arg
          break
        endif

        if( valnum(arg) != 0 & ymax = 'none' )
          ymax = arg
          break
        endif

        say 'syntax error: 'arg
        return

      endwhile
    endwhile
*  else
*** Usage 2 ***
* font size
*    fwidth  = subwrd( args, 5 )
*    fheight = subwrd( args, 6 )
*    fstep   = subwrd( args, 7 )
*    foffset = subwrd( args, 8 )
*  endif


*** get shade information ***
  if( cnum = -1 )
    'q shades'
    shdinfo = result
    if ( subwrd( shdinfo, 1 ) = 'None' ) 
      say 'Error in xcbar.gs: No shading information'
      return
    endif

* number of colors
    cnum = subwrd( shdinfo, 5 )
    if( cnum <= 0 )
      say 'Error in xcbar.gs: Number of color is zero'
      return
    endif

* color and (higher) levels
    i = 1
    while( i <= cnum )
      rec = sublin( shdinfo, i+1 )
      col.i = subwrd( rec, 1 )
      hi.i = subwrd( rec, 3 )
      i = i + 1
    endwhile

  endif


*** determine xmin, xmax, ymin, ymax if necessary ***
* ( following cbar.gs )
  if( xmin = 'none' | xmax = 'none' | ymin = 'none' | ymax = 'none' )
    'q gxinfo'
    sline = sublin( result, 2 )
    xsize = subwrd( sline, 4 )
    ysize = subwrd( sline, 6 )
    xline = sublin( result, 3 )
    yline = sublin( result, 4 )
    xlmin = subwrd( xline, 4 )
    xlmax = subwrd( xline, 6 )
    xlwid = xlmax - xlmin
    ylmin = subwrd( yline, 4 )
    ylmax = subwrd( yline, 6 )
    ylwid = ylmax - ylmin

*   vertical
    if( ylmin < 0.6 | xsize-xlmax > 1.5 )
      direction = 'vertical'
      xmin = xlmax + ( xsize - xlmax ) / 2 - 0.4
      xmax = xmin + 0.2

      y1wid = 0.5
      if ( y1wid * cnum > ysize * 0.8 ) 
        y1wid = ysize * 0.8 / cnum
      endif
      ymin = ysize / 2 - y1wid * cnum / 2
      ymax = ysize / 2 + y1wid * cnum / 2

*   horizontal
    else
      direction = 'horizontal'
      ymin = ylmin / 2
      ymax = ymin + 0.2

      x1wid = 0.8
      if ( x1wid * cnum > xsize * 0.8 )
        x1wid = xsize * 0.8 / cnum
      endif
      xmin = xsize / 2 - x1wid * cnum / 2
      xmax = xsize / 2 + x1wid * cnum / 2

    endif

  endif

*** determine direction if necessary ***
  if( direction != 'horizontal' & direction != 'vertical' )
    if( xmax - xmin >= ymax - ymin )
      direction = 'horizontal'
    else
      direction = 'vertical'
    endif
  endif
  if( direction = 'horizontal' )
    xdir = 1
    ydir = 0
  else
    xdir = 0
    ydir = 1
  endif

*** other default value *****
  if( edge != 'box' & edge != 'triangle' & edge != 'circle' )
    edge = 'box'
  endif

  if( valnum(fwidth) = 0 )
    fwidth = 0.12
  endif

  if( valnum(fheight) = 0 )
    fheight = 0.13
  endif

  if( valnum(fthickness) = 0 )
    fthickness = fheight * 40
  endif

  if( valnum(fstep) = 0 )
    fstep = 1
  endif

  if( valnum(foffset) = 0 )
    foffset = 0
*  foffset = fstep - 1
  endif

  if( line != 'on' & line != 'off' )
    line = 'off'
  endif



*** get other constant ***

* width of color bar
  xdif = xdir * (xmax-xmin) / cnum
  ydif = ydir * (ymax-ymin) / cnum


* position of the fonts
  xmoji = xmin + fwidth * ydir
  ymoji = ymin - fheight * xdir


*** draw ***
  i = 1
  x1 = xmin - xdif
  x2 = xmin * xdir + xmax * ydir
  y1 = ymin - ydif
  y2 = ymin * ydir + ymax * xdir

  while( i <= cnum )
    x1 = x1 + xdif
    x2 = x2 + xdif
    y1 = y1 + ydif
    y2 = y2 + ydif
    xmoji = x2 + ( 0.5 * fwidth ) * ydir
    ymoji = ( y1 - 0.5 * fheight ) * xdir + y2 * ydir

*    rec = sublin( shdinfo, i+1 )
*    col = subwrd( rec, 1 )
*    hi = subwrd( rec, 3 )
    'set line 'col.i
    'set strsiz 'fwidth' 'fheight

*** draw color bar ***
    if( edge = 'box' )
      'draw recf 'x1' 'y1' 'x2' 'y2
      if( line = 'on' ) ; drawrec( x1, y1, x2, y2 ) ; endif
    endif

    if( edge = 'triangle' )
      if( i != 1 & i != cnum )
        'draw recf 'x1' 'y1' 'x2' 'y2
        if( line = 'on' ) ; drawrec( x1, y1, x2, y2 ) ; endif
      endif

      if( direction = 'horizontal' )
        if( i = 1 )
          ymed = 0.5 * ( y1 + y2 )
          poly = x1' 'ymed' 'x2' 'y1' 'x2' 'y2
          'draw polyf 'poly
          if( line = 'on' ) ; drawpoly( poly ) ; endif
        endif
        if( i = cnum )
          ymed = 0.5 * ( y1 + y2 )
          poly = x1' 'y1' 'x1' 'y2' 'x2' 'ymed
          'draw polyf 'poly
          if( line = 'on' ) ; drawpoly( poly ) ; endif
        endif
      endif

      if( direction = 'vertical' )
        if( i = 1 )
          xmed = 0.5 * ( x1 + x2 )
          poly = xmed' 'y1' 'x1' 'y2' 'x2' 'y2
          'draw polyf 'poly
          if( line = 'on' ) ; drawpoly( poly ) ; endif
        endif
        if( i = cnum )
          xmed = 0.5 * ( x1 + x2 )
          poly = x1' 'y1' 'x2' 'y1' 'xmed' 'y2
          'draw polyf 'poly
          if( line = 'on' ) ; drawpoly( poly ) ; endif
        endif
      endif
    endif
*   end of triangle


    if( edge = 'circle' )
      if( i != 1 & i != cnum )
        'draw recf 'x1' 'y1' 'x2' 'y2
        if( line = 'on' ) ; drawrec( x1, y1, x2, y2 ) ; endif
      endif

      if( direction = 'horizontal' )
        if( i = 1 )
          xc = x1 + ( y2 - y1 )
          yc = y1 + 0.5 * ( y2 - y1 )
          radius = 0.5 * ( y2 - y1 )
          'draw recf 'xc' 'y1' 'x2' 'y2
          circle = circle( xc, yc, radius, 90, 270, 6 )
          'draw polyf 'circle
          if( line = 'on' ) ; drawpoly( x2' 'y2' 'circle' 'x2' 'y1 ) ; endif
        endif
        if( i = cnum )
          xc = x2 - ( y2 - y1 )
          yc = y1 + 0.5 * ( y2 - y1 )
          radius = 0.5 * ( y2 - y1 )
          'draw recf 'x1' 'y1' 'xc' 'y2
          circle = circle( xc, yc, radius, 270, 450, 6 )
          'draw polyf 'circle
          if( line = 'on' ) ; drawpoly(x1' 'y1' 'circle' 'x1' 'y2) ; endif
        endif
      endif

      if( direction = 'vertical' )
        if( i = 1 )
          xc = x1 + 0.5 * ( x2 - x1 )
          yc = y1 + 0.5 * ( x2 - x1 )
          radius = 0.5 * ( x2 - x1 )
          'draw recf 'x1' 'yc' 'x2' 'y2
          circle = circle( xc, yc, radius, 180, 360, 6 )
          'draw polyf 'circle
          if( line = 'on' ) ; drawpoly( x1' 'y2' 'circle' 'x2' 'y2 ) ; endif
        endif
        if( i = cnum )
          xc = x1 + 0.5 * ( x2 - x1 )
          yc = y2 - 0.5 * ( x2 - x1 )
          radius = 0.5 * ( x2 - x1 )
          'draw recf 'x1' 'y1' 'x2' 'yc
          circle = circle( xc, yc, radius, 0, 180, 6 )
          'draw polyf 'circle
          if( line = 'on' ) ; drawpoly( x2' 'y1' 'circle' 'x1' 'y1 ) ; endif
        endif
      endif
    endif
*   end of circle


* draw labels
*  if( i != cnum & i-math_int(i/fstep)*fstep = foffset)
    if( i != cnum & i-foffset > 0 & math_int((i-1-foffset)/fstep)*fstep = i-1-foffset  )

      if( direction = 'horizontal' )
        'set string 1 tc 'fthickness' 0'
      else
       'set string 1 l 'fthickness' 0'
      endif

      'draw string 'xmoji' 'ymoji' 'hi.i

    endif

    i = i + 1
  endwhile

return



*
* angle = 0  : x(+) direction
* angle = 90 : y(+) direction
*
function circle( xc, yc, radius, amin, amax, astep )
  circle=''

  angle = amin
  while( angle <= amax )
    x = xc + radius * math_cos( angle * 3.14 / 180.0 )
    y = yc + radius * math_sin( angle * 3.14 / 180.0 )
    circle = circle % x % ' ' % y % ' '
    angle = angle + astep
  endwhile

return ( circle )



function drawrec( xmin, ymin, xmax, ymax )
  drawpoly( xmin' 'ymin' 'xmin' 'ymax' 'xmax' 'ymax' 'xmax' 'ymin )
return


function drawpoly( args )
  xstart = subwrd( args, 1 )
  ystart = subwrd( args, 2 )

  xmin = xstart
  ymin = ystart

  i = 3
  while( 1 = 1 )
    xmax = subwrd( args, i )
    ymax = subwrd( args, i+1 )
    if( xmax = "" | ymax = "" ); break; endif

    'set cthick 1'
    'set line 1'
    'draw line 'xmin' 'ymin' 'xmax' 'ymax

    xmin = xmax
    ymin = ymax

    i = i + 2
  endwhile


  'set cthick 1'
  'set line 1'
  'draw line 'xmin' 'ymin' 'xstart' 'ystart

return


*
* help
*
function help()
  say ' Name:'
  say '   xcbar '_version' - Draw color bar at any position and size'
  say ' '
  say ' Usage:'
  say '   xcbar [ xmin xmax ymin ymax ]'
  say '         [ ( -fwidth | -fw ) value ]'
  say '         [ ( -fheight | -fh ) value ]'
  say '         [ ( -fthickness | -ft ) value ]'
  say '         [ ( -fstep | -fs ) value ]'
  say '         [ ( -foffset | -fo ) value ]'
  say '         [ ( -direction | -dir ) ( horizontal | h | vertical | v ) ]'
  say '         [ -edge ( box | triangle | circle ) ]'
  say '         [ -line ( on | off) ]'
  say '         [ -levcol c(1) l(1) c(2) level(2) ... l(cnum-1) c(cnum) ]'
  say ''
  say '     xmin       : color bar position (left side)'
  say '     xmax       : color bar position (right side)'
  say '     ymin       : color bar position (bottom side)'
  say '     ymax       : color bar position (top side)'
  say '                  without xmin, xmax, ymin or ymax, '
  say '                  position will be determined automatically'
  say '                  following cbar.gs manner.'
  say '     fwidth     : font width (default=0.12)'
  say '     fw         : same as fwidth'
  say '     fheight    : font height (default=0.13)'
  say '     fh         : same as fheight'
  say '     fthickness : font thickness (default=fheight*40)'
  say '     ft         : same as fthickness'
  say '     fstep      : label step (default=1)'
  say '     fs         : same as fstep'
  say '     foffset    : label offset for fstep (default=0)'
  say '     fo         : same as foffset'
  say '     direction  : horizontal ("h" in short) or vertical ("v" in short)'
  say '                  color bar (default=horizontal)'
  say '     dir        : same as direction'
  say '     edge       : shape of edge (default=box)'
  say '     line       : lines between each color box (default=off)'
  say '     c(1) l(1) c(2) level(2) ... l(cnum-1) c(cnum)'
  say '                : color numbers and levels. By using this option,'
  say '                  you can draw color bar without drawing figure'
  say ''
  say ' Note:'
  say '   [arg-name]       : specify if needed'
  say '   (arg1 | arg2)    : arg1 or arg2 must be specified'
  say ''
  say '    xcbar is based on cbar.gs'
  say ''
  say ' Copyright (C) 2009 Chihiro Kodama'
  say ' Distributed under GNU GPL (http://www.gnu.org/licenses/gpl.html)'
  say ''
return
