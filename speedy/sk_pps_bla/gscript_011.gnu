 set terminal tgif
 set terminal postscript eps enhanced color solid "Helvetica" 32
#  set terminal postscript eps "Helvetica" 32
 set size 1.3
# set size ratio 0.8 1.5
 set size ratio 0.7 2.0

 #set grid
 set key left top box lc rgb "black" lw 2


set size ratio 0.5 2.0
###########################################################################################################################################################
# (1)  :: YYYYMMDDHH
# (2)  :: STEPS 
# (3)  :: TIMER = STEPS*0.25
# (4)  :: EXPER RMSE
# (5)  :: EXPER SPRD
# (6)  :: EXPER RMSE 7-day mean (as of 20200415)
# (7)  :: EXPER SPRD 7-day mean (as of 20200415)
# (8)  :: EXPER Neff
# (9)  :: EXPER Neff 7-day mean (as of 20200415)

###########################################################################################################################################################
set ylabel "RMSE, Spread [K]"   font "Helvetica,42" offset 0.5,0
set xlabel "month/date in 1982" font "Helvetica,42" offset 0.5,0
###########################################################################################################################################################
set key right top
set title  "(a) [FG] RMSE (solid) & Spread (dashed) ; Temperature [K] at 500 hPa" font "Helvetica,40"
set yrange [0.05:0.55] ; set ytics 0, 0.10
set xrange [31.25:151.0] 
set xtics ('02/01'31.25, '03/01'59, '04/01'90, '05/01'120, '06/01'151, '07/01'181, '08/01'212, '09/01'243 ) font "Helvetica,33"
set output "./fig_011/time_reg2_rmse_lpf.eps"
plot \
  'gnu_000/time_reg2_2LAPF__M000040L0500IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:5 w l lw 10 lt 0 lc rgb "red"        title "",\
  'gnu_000/time_reg2_2LAPF__M000040L0400IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:5 w l lw 10 lt 0 lc rgb "blue"       title "",\
  'gnu_000/time_reg2_2LAPF__M000040L0300IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:5 w l lw 10 lt 0 lc rgb "dark-green" title "",\
  'gnu_000/time_reg2_LETKF_M000040L0500IADP.txt'                                                u 3:5 w l lw 10 lt 0 lc rgb "black"      title "",\
  'gnu_000/time_reg2_LETKF_M000040L0500IADP.txt'                                                u 3:4 w l lw 10 lt 1 lc rgb "black"      title "LETKF(h:500,adp-infl)",\
  'gnu_000/time_reg2_2LAPF__M000040L0500IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:4 w l lw 10 lt 1 lc rgb "red"        title "LPF(h:500,RTPS=1.00)",\
  'gnu_000/time_reg2_2LAPF__M000040L0400IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:4 w l lw 10 lt 1 lc rgb "blue"       title "LPF(h:400,RTPS=1.00)",\
  'gnu_000/time_reg2_2LAPF__M000040L0300IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:4 w l lw 10 lt 1 lc rgb "dark-green" title "LPF(h:300,RTPS=1.00)"
###
set output "./fig_011/time_raob_rmse_lpf.eps"
set yrange [0.05:1.05] ; set ytics 0, 0.10
plot \
  'gnu_000/time_raob_2LAPF__M000040L0400IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:5 w l lw 10 lt 0 lc rgb "dark-green" title "",\
  'gnu_000/time_raob_2LAPF__M000040L0600IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:5 w l lw 10 lt 0 lc rgb "blue"       title "",\
  'gnu_000/time_raob_2LAPF__M000040L0800IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:5 w l lw 10 lt 0 lc rgb "red"        title "",\
  'gnu_000/time_raob_LETKF_M000040L0500IADP.txt'                                                u 3:5 w l lw 10 lt 0 lc rgb "black"      title "",\
  'gnu_000/time_raob_LETKF_M000040L0500IADP.txt'                                                u 3:4 w l lw 10 lt 1 lc rgb "black"      title "LETKF(h:800,adp-infl)",\
  'gnu_000/time_raob_2LAPF__M000040L0400IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:4 w l lw 10 lt 1 lc rgb "dark-green" title "LPF(h:400,RTPS=1.00)",\
  'gnu_000/time_raob_2LAPF__M000040L0600IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:4 w l lw 10 lt 1 lc rgb "blue"       title "LPF(h:600,RTPS=1.00)",\
  'gnu_000/time_raob_2LAPF__M000040L0800IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:4 w l lw 10 lt 1 lc rgb "red"        title "LPF(h:800,RTPS=1.00)"

##############################
set yrange [0.11:0.26] ; set ytics 0, 0.04
set output "./fig_011/time_reg2_rmse_lpfgm_fgt0.00.eps"
plot \
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:5 w l lw 10 lt 0 lc rgb "magenta"    title "",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.80_woWSM_infKK_fac2.5_fgt0.00_rsmp010.txt'   u 3:5 w l lw 10 lt 0 lc rgb "red"        title "",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.80_woWSM_infKK_fac2.5_fgt0.00_rsmp005.txt'   u 3:5 w l lw 10 lt 0 lc rgb "blue"       title "",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.80_woWSM_infKK_fac2.5_fgt0.00_rsmp002.txt'   u 3:5 w l lw 10 lt 0 lc rgb "dark-green" title "",\
  'gnu_000/time_reg2_LETKF_M000040L0500IADP.txt'                                                u 3:5 w l lw 10 lt 0 lc rgb "black"      title "",\
  'gnu_000/time_reg2_LETKF_M000040L0500IADP.txt'                                                u 3:4 w l lw 10 lt 1 lc rgb "black"      title "LETKF(h:500,adp-infl)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:4 w l lw 10 lt 1 lc rgb "magenta"    title "LPFGM(h:500,RTPS=1.00,t=0.0,N_0=40)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.80_woWSM_infKK_fac2.5_fgt0.00_rsmp010.txt'   u 3:4 w l lw 10 lt 1 lc rgb "red"        title "LPFGM(h:500,RTPS=0.80,t=0.0,N_0=10)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.80_woWSM_infKK_fac2.5_fgt0.00_rsmp005.txt'   u 3:4 w l lw 10 lt 1 lc rgb "blue"       title "LPFGM(h:500,RTPS=0.80,t=0.0,N_0=5)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.80_woWSM_infKK_fac2.5_fgt0.00_rsmp002.txt'   u 3:4 w l lw 10 lt 1 lc rgb "dark-green" title "LPFGM(h:500,RTPS=0.80,t=0.0,N_0=2)"
#####
set key left top
set output "./fig_011/time_reg2_rmse_lpfgm_fgt1.00.eps"
plot \
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:5 w l lw 10 lt 0 lc rgb "magenta"    title "",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.50_woWSM_infKK_fac2.5_fgt1.00_rsmp010.txt'   u 3:5 w l lw 10 lt 0 lc rgb "red"        title "",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.50_woWSM_infKK_fac2.5_fgt1.00_rsmp005.txt'   u 3:5 w l lw 10 lt 0 lc rgb "blue"       title "",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.50_woWSM_infKK_fac2.5_fgt1.00_rsmp002.txt'   u 3:5 w l lw 10 lt 0 lc rgb "dark-green" title "",\
  'gnu_000/time_reg2_LETKF_M000040L0500IADP.txt'                                                u 3:5 w l lw 10 lt 0 lc rgb "black"      title "",\
  'gnu_000/time_reg2_LETKF_M000040L0500IADP.txt'                                                u 3:4 w l lw 10 lt 1 lc rgb "black"      title "LETKF(h:500,adp-infl)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:4 w l lw 10 lt 1 lc rgb "magenta"    title "LPFGM(h:500,RTPS=1.00,N_0=40)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.50_woWSM_infKK_fac2.5_fgt1.00_rsmp010.txt'   u 3:4 w l lw 10 lt 1 lc rgb "red"        title "LPFGM(h:500,RTPS=0.50,t=1.0,N_0=10)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.50_woWSM_infKK_fac2.5_fgt1.00_rsmp005.txt'   u 3:4 w l lw 10 lt 1 lc rgb "blue"       title "LPFGM(h:500,RTPS=0.50,t=1.0,N_0=5)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.50_woWSM_infKK_fac2.5_fgt1.00_rsmp002.txt'   u 3:4 w l lw 10 lt 1 lc rgb "dark-green" title "LPFGM(h:500,RTPS=0.50,t=1.0,N_0=2)"


########################################################################################################################################################### raob
set yrange [0.10:1.00] ; set ytics 0, 0.20
set key right top
set output "./fig_011/time_raob_rmse_lpfgm_fgt1.00.eps"
plot \
  'gnu_000/time_raob_3LPFGM_M000040L0800IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:5 w l lw 10 lt 0 lc rgb "magenta"    title "",\
  'gnu_000/time_raob_3LPFGM_M000040L0800IADP_RTPS0.60_woWSM_infKK_fac2.5_fgt1.00_rsmp010.txt'   u 3:5 w l lw 10 lt 0 lc rgb "red"        title "",\
  'gnu_000/time_raob_3LPFGM_M000040L0800IADP_RTPS0.60_woWSM_infKK_fac2.5_fgt1.00_rsmp005.txt'   u 3:5 w l lw 10 lt 0 lc rgb "blue"       title "",\
  'gnu_000/time_raob_3LPFGM_M000040L0800IADP_RTPS0.60_woWSM_infKK_fac2.5_fgt1.00_rsmp002.txt'   u 3:5 w l lw 10 lt 0 lc rgb "dark-green" title "",\
  'gnu_000/time_raob_LETKF_M000040L0800IADP.txt'                                                u 3:5 w l lw 10 lt 0 lc rgb "black"      title "",\
  'gnu_000/time_raob_LETKF_M000040L0800IADP.txt'                                                u 3:4 w l lw 10 lt 1 lc rgb "black"      title "LETKF(h:800,adp-infl)",\
  'gnu_000/time_raob_3LPFGM_M000040L0800IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:4 w l lw 10 lt 1 lc rgb "magenta"    title "LPFGM(h:800,RTPS=1.00,N_0=40)",\
  'gnu_000/time_raob_3LPFGM_M000040L0800IADP_RTPS0.60_woWSM_infKK_fac2.5_fgt1.00_rsmp010.txt'   u 3:4 w l lw 10 lt 1 lc rgb "red"        title "LPFGM(h:800,RTPS=0.60,t=1.0,N_0=10)",\
  'gnu_000/time_raob_3LPFGM_M000040L0800IADP_RTPS0.60_woWSM_infKK_fac2.5_fgt1.00_rsmp005.txt'   u 3:4 w l lw 10 lt 1 lc rgb "blue"       title "LPFGM(h:800,RTPS=0.60,t=1.0,N_0=5)",\
  'gnu_000/time_raob_3LPFGM_M000040L0800IADP_RTPS0.60_woWSM_infKK_fac2.5_fgt1.00_rsmp002.txt'   u 3:4 w l lw 10 lt 1 lc rgb "dark-green" title "LPFGM(h:800,RTPS=0.60,t=1.0,N_0=2)"

#  'gnu_000/time_raob_LETKF_M000040L0800IADP.txt'                                                u 3:5 w l lw 10 lt 0 lc rgb "black"      title "",\
#  'gnu_000/time_raob_LETKF_M000040L0800IADP.txt'                                                u 3:4 w l lw 10 lt 1 lc rgb "black"      title "LETKF(h:800,adp-infl)",\


###########################################################################################################################################################
set ylabel "Neff"   font "Helvetica,42" offset 0.5,0
set title  "(b) [FG] Effective Particle Size :: Temperature [K] at 500 hPa" font "Helvetica,40"
###########################################################################################################################################################
set key center top
set yrange [20:40] ; set ytics 0, 2
set output "./fig_011/time_reg2_peff_lpf.eps"
plot \
  'gnu_000/time_reg2_2LAPF__M000040L0300IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:8 w l lw 10 lt 1 lc rgb "dark-green" title "LPF(h:300,RTPS=1.00)",\
  'gnu_000/time_reg2_2LAPF__M000040L0400IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:8 w l lw 10 lt 1 lc rgb "blue"       title "LPF(h:400,RTPS=1.00)",\
  'gnu_000/time_reg2_2LAPF__M000040L0500IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:8 w l lw 10 lt 1 lc rgb "red"        title "LPF(h:500,RTPS=1.00)"
##############################
set yrange [30:37] ; set ytics 0, 2
set output "./fig_011/time_reg2_peff_lpfgm_fgt1.00.eps"
plot \
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:8 w l lw 10 lt 1 lc rgb "magenta"    title "LPFGM(h:500,RTPS=1.00,N_0=40)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.50_woWSM_infKK_fac2.5_fgt1.00_rsmp010.txt'   u 3:8 w l lw 10 lt 1 lc rgb "red"        title "LPFGM(h:500,RTPS=0.50,t=1.0,N_0=10)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.50_woWSM_infKK_fac2.5_fgt1.00_rsmp005.txt'   u 3:8 w l lw 10 lt 1 lc rgb "blue"       title "LPFGM(h:500,RTPS=0.50,t=1.0,N_0=5)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.50_woWSM_infKK_fac2.5_fgt1.00_rsmp002.txt'   u 3:8 w l lw 10 lt 1 lc rgb "dark-green" title "LPFGM(h:500,RTPS=0.50,t=1.0,N_0=2)"

##############################
set yrange [30:37] ; set ytics 0, 2
set output "./fig_011/time_reg2_peff_lpfgm_fgt0.00.eps"
plot \
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:8 w l lw 10 lt 1 lc rgb "magenta"    title "LPFGM(h:500,RTPS=1.00,t=0.0,N_0=40)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.80_woWSM_infKK_fac2.5_fgt0.00_rsmp010.txt'   u 3:8 w l lw 10 lt 1 lc rgb "red"        title "LPFGM(h:500,RTPS=0.80,t=0.0,N_0=10)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.80_woWSM_infKK_fac2.5_fgt0.00_rsmp005.txt'   u 3:8 w l lw 10 lt 1 lc rgb "blue"       title "LPFGM(h:500,RTPS=0.80,t=0.0,N_0=5)",\
  'gnu_000/time_reg2_3LPFGM_M000040L0500IADP_RTPS0.80_woWSM_infKK_fac2.5_fgt0.00_rsmp002.txt'   u 3:8 w l lw 10 lt 1 lc rgb "dark-green" title "LPFGM(h:500,RTPS=0.80,t=0.0,N_0=2)"

##########################################################
set key left bottom
set yrange [18:32] ; set ytics 0, 2
set output "./fig_011/time_raob_peff_lpfgm_fgt1.00.eps"
plot \
  'gnu_000/time_raob_3LPFGM_M000040L0800IADP_RTPS1.00_woWSM_infKK_fac2.5_fgt0.00_rsmp040.txt'   u 3:8 w l lw 10 lt 1 lc rgb "magenta"    title "LPFGM(h:800,RTPS=1.00,N_0=40)",\
  'gnu_000/time_raob_3LPFGM_M000040L0800IADP_RTPS0.60_woWSM_infKK_fac2.5_fgt1.00_rsmp010.txt'   u 3:8 w l lw 10 lt 1 lc rgb "red"        title "LPFGM(h:800,RTPS=0.60,t=1.0,N_0=10)",\
  'gnu_000/time_raob_3LPFGM_M000040L0800IADP_RTPS0.60_woWSM_infKK_fac2.5_fgt1.00_rsmp005.txt'   u 3:8 w l lw 10 lt 1 lc rgb "blue"       title "LPFGM(h:800,RTPS=0.60,t=1.0,N_0=5)",\
  'gnu_000/time_raob_3LPFGM_M000040L0800IADP_RTPS0.60_woWSM_infKK_fac2.5_fgt1.00_rsmp002.txt'   u 3:8 w l lw 10 lt 1 lc rgb "dark-green" title "LPFGM(h:800,RTPS=0.60,t=1.0,N_0=2)"