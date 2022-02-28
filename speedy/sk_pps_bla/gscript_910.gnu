 set terminal tgif
 set terminal postscript eps enhanced color solid "Helvetica" 32
#  set terminal postscript eps "Helvetica" 32
 set size 1.3
# set size ratio 0.8 1.5
 set size ratio 0.7 2.0

# set grid
# set key left top box lc rgb "black" lw 2


set size ratio 0.5 2.0

###########################################################################################################################################################
set ylabel "Absolute Error"                     font "Helvetica,40" offset 0.5,0
set xlabel "the Number of Monte-Carlo Samples"  font "Helvetica,40" offset 0.5,0
set key right top
###########################################################################################################################################################
set title  "(a) Abs. Error b/w Monte-Carlo-based T and prescribed w (M=040)" font "Helvetica,40"
#set yrange [0.37:0.62] ; set ytics 0, 0.05
set xrange [0.0:500.0] 
set output "./fig_910/abs_error_vec.eps"
plot \
  'out_900/sampled-error_M000040NTIME00010000NCASE00001000.txt'   u 2:7:8 w filledcurves lc rgb "orange" title "min-max",\
  'out_900/sampled-error_M000040NTIME00010000NCASE00001000.txt'   u 2:3 w l lw 10 lt 1   lc rgb "red"    title "average"
################################################################################################################
set title  "(b) Abs. Error b/w Monte-Carlo-based T and T with 10000 samples (M=040)" font "Helvetica,40"
#set yrange [0.37:0.62] ; set ytics 0, 0.05
set xrange [0.0:500.0] 
set output "./fig_910/abs_error_mtx.eps"
plot \
  'out_900/sampled-error_M000040NTIME00010000NCASE00001000.txt'   u 2:13:14 w filledcurves lc rgb "orange" title "min-max",\
  'out_900/sampled-error_M000040NTIME00010000NCASE00001000.txt'   u 2:9     w l lw 10 lt 1 lc rgb "red"    title "average"

#  'out_900/sampled-error_M000040NTIME00010000NCASE00001000.txt'   u 2:5:6 w filledcurves lc rgb "orange" title "",\
#  
