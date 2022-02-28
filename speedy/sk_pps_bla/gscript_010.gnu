 set terminal tgif
 set terminal postscript eps enhanced color solid "Helvetica" 32
#  set terminal postscript eps "Helvetica" 32
 set size 1.3
# set size ratio 0.8 1.5
 set size ratio 0.7 2.0

 #set grid
 set key left top box lc rgb "black" lw 2

###########################################################################################################################################################
# set key box
# set key below
 set size ratio 0.7 2.0

set xrange [-0.05:1.05] ; set xtics 0, 0.10
#set xtics ('GL'5, 'NH'15.0, 'TR'25.0, 'SH'35.0 ) font "Helvetica,33"
###################################################################
set ylabel "RMSE [K]"                           font "Helvetica,42" offset 0.5,0
set xlabel "alpha (RTPS parameter )"            font "Helvetica,42" offset 0.5,0
###################################################################
set xrange [0.75:1.25]
set yrange [0.00:5.00] ; set ytics 0, 0.50
set title  "(a) SPEEDY-LPF ; RMSE (FG) Temperature [K] at 500 hPa" font "Helvetica,42"
set key right top
set output "./fig_010/rtps_reg2_2LAPF__all_Hxxxx.eps"
plot \
  0.1450 w l  lw 10  lt 0         lc rgb "black" title "CTRL(LETKF;h500;adp-infl)",\
  'gnu_cat/tave_reg2_2LAPF__M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "red"        title "LPF(h:500)",\
  'gnu_cat/tave_reg2_2LAPF__M000040L0400IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "blue"       title "LPF(h:400)",\
  'gnu_cat/tave_reg2_2LAPF__M000040L0300IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "dark-green" title "LPF(h:300)"
###
set yrange [0.00:0.50] ; set ytics 0, 0.10
set ylabel""
set title  "(b) SPEEDY-LPF ; RMSE (FG) Temperature [K] at 500 hPa" font "Helvetica,42"
set key left bottom
set output "./fig_010/rtps_reg2_2LAPF__lim_Hxxxx.eps"
plot \
  0.1450 w l  lw 10  lt 0         lc rgb "black" title "CTRL(LETKF)",\
  'gnu_cat/tave_reg2_2LAPF__M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "red"        title "LPF(h:500)",\
  'gnu_cat/tave_reg2_2LAPF__M000040L0400IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "blue"       title "LPF(h:400)",\
  'gnu_cat/tave_reg2_2LAPF__M000040L0300IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "dark-green" title "LPF(h:300)"

set ylabel "RMSE [K]"                           font "Helvetica,42" offset 0.5,0
##########
set yrange [0.00:5.00] ; set ytics 0, 0.40
set title  "(a) SPEEDY ; RMSE (FG) Temperature [K] at 500 hPa" font "Helvetica,42"
set output "./fig_010/rtps_raob_2LAPF__all_Hxxxx.eps"
plot \
  0.1450 w l  lw 10  lt 0         lc rgb "black" title "CTRL(LETKF;h500;adp-infl)",\
  'gnu_cat/tave_raob_2LAPF__M000040L0800IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "red"        title "LPF(h:800)",\
  'gnu_cat/tave_raob_2LAPF__M000040L0600IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "blue"       title "LPF(h:600)",\
  'gnu_cat/tave_raob_2LAPF__M000040L0400IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "dark-green" title "LPF(h:400)"
###
set yrange [0.00:0.50] ; set ytics 0, 0.10
set title  "(b) SPEEDY ; RMSE (FG) Temperature [K] at 500 hPa" font "Helvetica,42"
set output "./fig_010/rtps_raob_2LAPF__lim_Hxxxx.eps"
plot \
  0.1450 w l  lw 10  lt 0         lc rgb "black" title "CTRL(LETKF)",\
  'gnu_cat/tave_raob_2LAPF__M000040L0800IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "red"        title "LPF(h:800)",\
  'gnu_cat/tave_raob_2LAPF__M000040L0600IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "blue"       title "LPF(h:600)",\
  'gnu_cat/tave_raob_2LAPF__M000040L0400IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "dark-green" title "LPF(h:400)"


###################################################################
set xrange [0.00:1.00]
set yrange [0.00:0.40] ; set ytics 0, 0.10
set title  "(b) SPEEDY-LPFGM ; RMSE (FG) Temperature [K] at 500 hPa" font "Helvetica,42"
set output "./fig_010/rtps_reg2_3LPFGM_tau1.0_lim.eps"
plot \
  0.1450 w l  lw 10  lt 0         lc rgb "black" title "CTRL(LETKF;h500;adp-infl)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp010_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 0  lw 3 pt 7 ps 2 lc rgb "red"        title "LPFGM(h:500,N_0=10,t=1.0)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp005_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 0  lw 3 pt 7 ps 2 lc rgb "blue"       title "LPFGM(h:500,N_0=5,t=1.0)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp002_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 0  lw 3 pt 7 ps 2 lc rgb "dark-green" title "LPFGM(h:500,N_0=2,t=1.0)"

###
set title  "(a) SPEEDY-LPFGM ; RMSE (FG) Temperature [K] at 500 hPa" font "Helvetica,42"
set output "./fig_010/rtps_reg2_3LPFGM_tau0.0_lim.eps"
plot \
  0.1450 w l  lw 10  lt 0         lc rgb "black" title "CTRL(LETKF;h500;adp-infl)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "magenta"    title "LPFGM(h:500,N_0=40)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp010_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "red"        title "LPFGM(h:500,N_0=10,t=0.0)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp005_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "blue"       title "LPFGM(h:500,N_0=5,t=0.0)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp002_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "dark-green" title "LPFGM(h:500,N_0=2,t=0.0)"


###################################################################
set key right top
set xrange [0.00:1.00]
set yrange [0.40:5.00] ; set ytics 0, 0.50
set title  "(a) SPEEDY-LPFGM ; RMSE (FG) Temperature [K] at 500 hPa" font "Helvetica,42"
set output "./fig_010/rtps_raob_3LPFGM_tau1.0_all.eps"
plot \
  0.4762 w l  lw 10  lt 0         lc rgb "black" title "CTRL(LETKF;h800;adp-infl)",\
  'gnu_cat/tave_raob_3LPFGM_M000040L0800IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 7 ps 2 lc rgb "magenta"        title "LPFGM(h:800,N_0=40)",\
  'gnu_cat/tave_raob_3LPFGM_M000040L0800IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp010_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 0  lw 3 pt 7 ps 2 lc rgb "red"        title "LPFGM(h:800,N_0=10,t=0.0)",\
  'gnu_cat/tave_raob_3LPFGM_M000040L0800IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp005_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 0  lw 3 pt 7 ps 2 lc rgb "blue"       title "LPFGM(h:800,N_0=5,t=0.0)",\
  'gnu_cat/tave_raob_3LPFGM_M000040L0800IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp002_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 0  lw 3 pt 7 ps 2 lc rgb "dark-green" title "LPFGM(h:800,N_0=2,t=0.0)"
#
set yrange [0.40:0.80] ; set ytics 0, 0.10
set title  "(b) SPEEDY-LPFGM ; RMSE (FG) Temperature [K] at 500 hPa" font "Helvetica,42"
set output "./fig_010/rtps_raob_3LPFGM_tau1.0_lim.eps"
plot \
  0.4762 w l  lw 10  lt 0         lc rgb "black" title "CTRL(LETKF;h800;adp-infl)",\
  'gnu_cat/tave_raob_3LPFGM_M000040L0800IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp010_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 0  lw 3 pt 7 ps 2 lc rgb "red"        title "LPFGM(h:800,N_0=10,t=0.0)",\
  'gnu_cat/tave_raob_3LPFGM_M000040L0800IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp005_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 0  lw 3 pt 7 ps 2 lc rgb "blue"       title "LPFGM(h:800,N_0=5,t=0.0)",\
  'gnu_cat/tave_raob_3LPFGM_M000040L0800IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp002_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 0  lw 3 pt 7 ps 2 lc rgb "dark-green" title "LPFGM(h:800,N_0=2,t=0.0)"


quit
stop

##########
set output "./fig_010/rtps_reg2_3LPFGM_tau_all.eps"
set yrange [0.00:5.00] ; set ytics 0, 0.40
plot \
  0.1450 w l  lw 10  lt 0         lc rgb "black" title "CTRL(LETKF;h500;adp-infl)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp005_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 0  lw 3 pt 2 ps 2 lc rgb "red"        title "LPFGM(h:500,t=0,N_0=5)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp010_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 0  lw 3 pt 2 ps 2 lc rgb "blue"       title "LPFGM(h:500,t=0,N_0=10)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp020_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 0  lw 3 pt 2 ps 2 lc rgb "dark-green" title "LPFGM(h:500,t=0,N_0=20)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp005_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "red"        title "LPFGM(h:500,t=1,N_0=5)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp010_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "blue"       title "LPFGM(h:500,t=1,N_0=10)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp020_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "dark-green" title "LPFGM(h:500,t=1,N_0=20)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt1.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "magenta"    title "LPFGM(h:500,t=1,N_0=40)"
##########
set output "./fig_010/rtps_reg2_3LPFGM_fac_all.eps"
set xrange [0.75:1.25]
set yrange [0.00:5.00] ; set ytics 0, 0.40
set yrange [0.00:0.50] ; set ytics 0, 0.10
plot \
  0.1450 w l  lw 10  lt 0         lc rgb "black" title "CTRL(LETKF;h500;adp-infl)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac1.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "red"        title "LPFGM(h:500,fac=1.5,N_0=40)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.0_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "blue"       title "LPFGM(h:500,fac=2.0,N_0=40)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac2.5_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "dark-green" title "LPFGM(h:500,fac=2.5,N_0=40)",\
  'gnu_cat/tave_reg2_3LPFGM_M000040L0500IADP_RTPSx.xx_woWSM_infKK_fac3.0_fgt0.00_rsmp040_SMP1982030100-1982053118.txt'   u 2:4 w lp lt 1  lw 3 pt 2 ps 2 lc rgb "magenta"    title "LPFGM(h:500,fac=3.0,N_0=40)"









