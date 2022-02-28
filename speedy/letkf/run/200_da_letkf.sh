#!/bin/bash
#=======================================================================
# letkf_cycle.sh
#   To run the SPEEDY-LETKF cycle in parallel computing environment
#=======================================================================
set -e
CDIR=`pwd`
#----------------------INITIALIZE PARAMETERS-----------------------------#
export NODE_list='n01'
export NODE=10                        #Nodes for parallel computing
export MPIRUN="mpirun"
LUTDIR='/../../src_LTERPLUT/output/'
#========================================================================== (1) DA OPTIONS
SYS='SYS20200810_RELAXATION_'
MEMTBL='20 30'
OBSTBL='raob' # raob, reg2, raob
HSGTBL='800'
#========================================================================== (2) WI OPTIONS
SMPTBL='1'      # sampling interval 
LUTTBL='BLLUT'  # LUT of WI
STHTBL='0'      # 0: no weight smoothing (default), 1: weight smoothing
#========================================================================== (3) DATA ASSIMILATION TYPE
DASTBL='3'
  # 0 :: LETKF    ;                    mean update :: LETKF,   ptb update :: LETKF
  # 1 :: DEBUG    ;                    mean update :: LPF,     ptb update :: LETKF                 # just for debug
  # 2 :: LAPF     ;                    mean update :: LPF  ,   ptb update :: LAPF      resampling
  # 3 :: LPFGM    ; LETKF_particle  +  mean update :: N/A  ,   ptb update :: LPFGM w   resampling

RTPTBL='2'      # 0: no relaxation  ;  1: RTPP  ;  2:RTPS (default)
RLXTBL='0.60'   # relaxation parameter alpha (RTPP)
				# 0: no relaxation   1: perfect relaxation i.e., dXa=dXb

#>> CAUTION <<# no relaxation is recommended for LETKF

#========================================================================== (3) LPF OPTIONS (usually fixed)

LFACTOR=2.5     # (usually fixed) factor of background error covariance (only for dastype=3 or 4)
RSMTBL='2'      # 0: shared resampling noise (Roland et al. 2019; MWR)
				# 1: KK's Monte-Carlo Resampling (repeat M-times)
				# 2: Marginal Particle Filter (MPF)
				#    (This parameter is currently used only for DASTYPE 2 or 3)
FGTTBL='1.00'   # 0.00: no forget, 1.00: perfect forget (default for LPFGM of raob)
SIRPEF=2       # effective particle size for resampling (SIR)
#==========================================================================
SYYYY=1982 ; SMM=01 ; SDD=01 ; SHH=00  # start of experiments 
EYYYY=1982 ; EMM=05 ; EDD=31 ; EHH=18  # end   of ecperiments
#==========================================================================

############################## parameter needed to be editted (end) ##############################
DISK=0 # 0:HDD, 1:RAM
NNODE=1
PPN=$NODE
PE=$NODE
export MPI_XPMEM_ENABLED=disabled
#------------------------------------------------------------------------------------------------#
NOSAVEENS=1 #0:save all members, 1:save only mean/sprd 
date +"%Y/%m/%d %H:%M:%S"
set_START(){
        local dummy
        read START dummy < /proc/uptime
}

get_ELAPS() {
        local dummy
        read END dummy < /proc/uptime
        let ELAPS=${END/./}0-${START/./}0
        ELAPS=`echo "scale=8; ${ELAPS}*0.001" | bc`
}
################################################################################## exp loop start (0) :: local particle filter
for DASTYPE in ${DASTBL}; do
for STH in ${STHTBL} ; do
for FGT in ${FGTTBL} ; do
for RSM in ${RSMTBL} ; do
for RTP in ${RTPTBL} ; do
for RLX in ${RLXTBL} ; do
################################################################################## exp loop start (1) :: weight interpolation
for LUT in ${LUTTBL} ; do
for SMP in ${SMPTBL} ; do
################################################################################## exp loop start (2) :: general
for MEMBER in ${MEMTBL} ; do
for OBS    in ${OBSTBL} ; do
for HSG    in ${HSGTBL} ; do
  HSIG=$HSG

  M=$( printf %06i ${MEMBER} )
  H=$( printf %04i ${HSIG}   )
  #==========================================>>>>> weight interpolation  (str)
  CSMP=$( printf %02i ${SMP} )
  
  if [ "$SMP" = "1" ]; then
    #CURRENT_TEST='_CONTROL'
    CURRENT_TEST=''
    LOGIC_INT=".false."
  else 
    LOGIC_INT=".true."
    #===> LUT-BASED WEIGHT INTERPOLATION
    CURRENT_TEST='_'$LUT$CSMP
    if [ $LUT = "BLLUT" ] ; then LUTNAM=$LUTDIR/interp_bilin${CSMP}.bin ; fi
    if [ $LUT = "ADJST" ] ; then LUTNAM=$LUTDIR/interp_adjst${CSMP}.bin ; fi
    if [ $LUT = "ADREV" ] ; then LUTNAM=$LUTDIR/interp_adjst${CSMP}.bin ; fi
    if [ $LUT = "REGUL" ] ; then LUTNAM=$LUTDIR/interp_regul${CSMP}.bin ; fi
    if [ $LUT = "RGREV" ] ; then LUTNAM=$LUTDIR/interp_regul${CSMP}.bin ; fi
    if [ $LUT = "WCUNI" ] ; then LUTNAM=$LUTDIR/interp_wcuni${CSMP}_${OBS}H${H}.bin ; fi
    if [ $LUT = "WCREG" ] ; then LUTNAM=$LUTDIR/interp_wcreg${CSMP}_${OBS}H${H}.bin ; fi

    if [ ! -e $LUTNAM ]; then
      echo "  error, LUTFILE NOT FOUND :: " $LUTNAM
      exit
    fi  
  fi

  #===> relaxation
  if [ $RTP = "0" ] ; then CURRENT_TEST=${CURRENT_TEST}'_RTPx'${RLX} ; fi
  if [ $RTP = "1" ] ; then CURRENT_TEST=${CURRENT_TEST}'_RTPP'${RLX} ; fi
  if [ $RTP = "2" ] ; then CURRENT_TEST=${CURRENT_TEST}'_RTPS'${RLX} ; fi        

  #===> weigt smoother
  if [ "$STH" = "0" ]; then CURRENT_TEST=${CURRENT_TEST}'_woWSM' ; fi
  if [ "$STH" = "1" ]; then CURRENT_TEST=${CURRENT_TEST}'_wwWSM' ; fi


  #===> LPF 
  if [ $DASTYPE -ge 2 ]; then
    if [ "$RSM" = "0" ]; then CURRENT_TEST=${CURRENT_TEST}'_infRP' ; fi
    if [ "$RSM" = "1" ]; then CURRENT_TEST=${CURRENT_TEST}'_infKK' ; fi
    if [ "$RSM" = "2" ]; then CURRENT_TEST=${CURRENT_TEST}'_infMP' ; fi
    CURRENT_TEST=${CURRENT_TEST}'_fac'${LFACTOR} # gamma of LPFGM

    CURRENT_TEST=${CURRENT_TEST}'_fgt'${FGT}     # forgetting factor      
    SIRPEF=$( printf %03i ${SIRPEF} )
    CURRENT_TEST=${CURRENT_TEST}'_rsmp'${SIRPEF} # resampling parameter
  fi

  #===> weight treatment
  if [ "$STH" = "0" ]; then LOGIC_STH=".false." ; fi
  if [ "$STH" = "1" ]; then LOGIC_STH=".true."  ; fi

  #===> DAS NAME
                               DASNAME="NULL__"
  if [ $DASTYPE -eq 0 ] ; then DASNAME="0LETKF_" ; fi
  if [ $DASTYPE -eq 1 ] ; then DASNAME="1LPFck_" ; fi
  if [ $DASTYPE -eq 2 ] ; then DASNAME="2LAPF__" ; fi
  if [ $DASTYPE -eq 3 ] ; then DASNAME="3LPFGM_" ; fi

  #==========================================>>>>> fixed descriptions  (str)
  EXP=M${M}L${H}IADP${CURRENT_TEST}
  MEM_list=$( seq 1 $M )  
  EXP=$SYS$DASNAME$EXP


### directory settings
cd $CDIR
cd ../..
SPEEDY=`pwd`
OUTPUT=$SPEEDY/DATA/$OBS/$EXP # data directory
OBSDIR=$SPEEDY/DATA/$OBS/obs  # obs data directory
TMPDIR=$SPEEDY/DATA/tmp_${OBS}_${EXP}/letkf # work directory
if test $DISK -eq 0 ; then
  TMPDIR=$SPEEDY/DATA/tmp_${OBS}_${EXP}/letkf # work directory
elif test $DISK -eq 1 ; then
  TMPDIR=/dev/shm/kotsuki_${OBS}_${EXP}/letkf   # work directory
fi

LETKF=letkf020.m01
OBSOPE=obsope.s01
OBSOPE=obsope_MPI.s01 # MPI'ed
#if [ -d $OUTPUT ] ; then rm -r $OUTPUT ; fi
if [ -d $OUTPUT ] ; then 
  #DEBUG#rm -r $OUTPUT
  echo "error :: stop not to overwrite results" $OUTPUT
  exit
fi

#==> add by kotsuki
INTDIR=$SPEEDY/DATA/init_ensemble
#if [ $OBS = "raob" ] ; then INTDIR=$SPEEDY/DATA/raob/init_ensemble_M000040L0800IADP ; fi
#if [ $OBS = "reg2" ] ; then INTDIR=$SPEEDY/DATA/reg2/init_ensemble_M000040L0500IADP ; fi

###if [ ! -d $OUTPUT ] ; then cp -ar $INTDIR $OUTPUT ; fi
mkdir -p $OUTPUT/gues 
mkdir -p $OUTPUT/wvec # SK
cd $OUTPUT
  cp -ar $INTDIR/anal      ./
  cp -ar $INTDIR/anal_f    ./
  cp -ar $INTDIR/infl_mul  ./
  cp -ar $INTDIR/peff_lpf  ./
  cp -ar $INTDIR/rsmp_lpf  ./
  cp -ar $INTDIR/asis_lpf  ./
  cp -ar $INTDIR/log       ./
  cp -ar $INTDIR/gues/mean ./gues/
  cp -ar $INTDIR/gues/sprd ./gues/
  cp -ar $INTDIR/gues/kldv ./gues/
cd $OUTPUT/gues
  for MEM in $MEM_list ; do
    MEM=$( printf %06i $MEM )
    cp -ar   $INTDIR/gues/$MEM ./
  done

### directory settings
#CDIR=`pwd`
echo $EXP
if test $DISK -eq 0 ; then
  TIMELOG="$CDIR/stdout_${OBS}_${EXP}_NODE${NODE}_HDD"
  echo $EXP with HDD > $TIMELOG
elif test $DISK -eq 1 ; then
  TIMELOG="$CDIR/stdout_${OBS}_${EXP}_NODE${NODE}_RAM"
  echo $EXP with RAM_DISK > $TIMELOG
fi
CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "Job started at $CTIME" >> $TIMELOG


cd $CDIR
### NODE info files
if test -e nodefile/$EXP ; then rm -fr nodefile/$EXP ; fi
mkdir -p nodefile/$EXP
# for LETKF
NODE_info=nodefile/$EXP/nodefile_letkf
rm -f $NODE_info ; touch $NODE_info
for n in $NODE_list ; do
  p=1
  while test $p -le $PPN ; do
    echo $n >> $NODE_info
    p=`expr $p + 1`
  done
done

#-----------------------------------------------------------------------
# Usually do not modify below
#-----------------------------------------------------------------------
source $SPEEDY/../common/timeinc.sh
  IYYYY=$SYYYY  ;  IMM=$SMM  ;  IDD=$SDD  ;  IHH=$SHH
  SB=$SPEEDY/model/data/bc/t30/clim
  SC=$SPEEDY/model/data/bc/t30/anom
#
# Work directory
#
rm -rf $TMPDIR
mkdir -p $TMPDIR
  cd $TMPDIR
  cp $CDIR/110_ensfcst_MPI.sh           ./
  cp $CDIR/120_ensemble_mean_MPI.sh     ./
  cp $CDIR/../ensemble_mean_MPI.exe     ./ 

  cp $SPEEDY/letkf/$LETKF               ./
  cp $SPEEDY/letkf/$OBSOPE              ./
  cp $CDIR/nodefile/$EXP/nodefile_letkf ./

  # Build SPEEDY model
  echo '>>> BEGIN BUILDING SPEEDY MODEL'
  cp $SPEEDY/model/run/imp_MPI_rank2.exe ./imp.exe
  echo '>>>'


### inputs
set_START
if test $DISK -eq 0 ; then
  ln -s $OUTPUT/anal     anal
  ln -s $OUTPUT/anal_f   anal_f
  ln -s $OUTPUT/gues     gues
  ln -s $OUTPUT/infl_mul infl_mul
  ln -s $OUTPUT/peff_lpf peff_lpf
  ln -s $OUTPUT/rsmp_lpf rsmp_lpf
  ln -s $OUTPUT/log      log
  ln -s $OBSDIR obs
  cp $SPEEDY/common/orography_t30.dat orography_t30.dat
  cp $SB/orog_lsm_alb.t30.grd         fort.20
  cp $SB/sst_8190clim.t30.sea.grd     fort.21
  cp $SB/seaice_8190clim.t30.sea.grd  fort.22
  cp $SB/skt_8190clim.t30.land.grd    fort.23
  cp $SB/sndep_8190clim.t30.land.grd  fort.24
  cp $SB/veget.t30.land.grd           fort.25
  cp $SB/soilw_8190clim.t30.land.grd  fort.26
  cp $SC/sst_anom_7990.t30.grd        fort.30
  CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "copy ended: $CTIME" >> $TIMELOG

  #==> SK
  ln -s $LUTNAM          BILINLUT.bin
  ln -s $OUTPUT/wvec     wvec
  ln -s $OUTPUT/asis_lpf asis_lpf

elif test $DISK -eq 1 ; then
  echo ">>> Compress files for STGIN"
  if test -e tmp_stgin ; then rm -rf tmp_stgin ; fi
  mkdir tmp_stgin
  cd tmp_stgin

    # Compress initial
    for MEM in $MEM_list ; do
      tMEM=$((MEM - 1))
      MEM=$( printf %06i $MEM )
  #   check=`expr $tMEM / $PPN`
  #   NODE_NUM=`expr $check % $NNODE + 1`
      check=$((tMEM / PPN))
      NODE_NUM=$((check % NNODE + 1))
      n=`echo $NODE_list | cut -d ' ' -f $NODE_NUM`

      mkdir -p $n/anal/$MEM
      mkdir -p $n/anal_f/$MEM
      mkdir -p $n/gues/$MEM
      echo " >> MEMBER:$MEM for $n ($NODE_NUM)"
      cp $OUTPUT/gues/$MEM/$SYYYY$SMM$SDD$SHH.grd $n/gues/$MEM
    done
  for n in $NODE_list ; do
    cd $n
    tar cf ../gues_$n.tar anal anal_f gues
    cd ../
  done
  CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "  initial ended:     $CTIME" >> $TIMELOG

  # Compress inputs
  echo " >> INPUTs"
  cp $SPEEDY/common/orography_t30.dat orography_t30.dat
  cp $SB/orog_lsm_alb.t30.grd         fort.20
  cp $SB/sst_8190clim.t30.sea.grd     fort.21
  cp $SB/seaice_8190clim.t30.sea.grd  fort.22
  cp $SB/skt_8190clim.t30.land.grd    fort.23
  cp $SB/sndep_8190clim.t30.land.grd  fort.24
  cp $SB/veget.t30.land.grd           fort.25
  cp $SB/soilw_8190clim.t30.land.grd  fort.26
  cp $SC/sst_anom_7990.t30.grd        fort.30
  tar cf other.tar orography_t30.dat  fort.??
  pwd
  echo $LUTNAM
  if [ "$SMP" -ne "1" ]; then
    echo "copying???"
    cp $LUTNAM                        BILINLUT.bin
    tar rf other.tar                  BILINLUT.bin
  fi
  CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "  copy ended:        $CTIME" >> $TIMELOG

  # Compress obs
  echo " >> OBSERVATION"
  mkdir obs
  IYMD=$SYYYY$SMM$SDD
  while test $IYMD -le $EYYYY$EMM$EDD ; do
    for n in $NODE_list ; do
      cp $SPEEDY/DATA/$OBS/obs/$IYMD??.dat obs
    done
    IYMD=`date -d "$IYMD 1 day" +%Y%m%d`
  done
  tar rf other.tar obs
  CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "  observation ended:$CTIME" >> $TIMELOG

  # Compress exe
  echo " >> EXE_files"
  cp ../*.exe          .
  cp ../obsope_MPI.s01 .
  cp ../nodefile_letkf .
  tar rf other.tar *.exe obsope_MPI.s01 nodefile_letkf
  CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "  exe files ended:  $CTIME" >> $TIMELOG

  # STAGE IN
  echo '>>>'
  echo '>>> STAGE IN'
  IYMD=$SYYYY$SMM$SDD
  PYMDH=`date -d "$IYMD $IHH:00 6 hour ago" +%Y%m%d%H`
  NODE_NUM=1
  for n in $NODE_list ; do
    n=`echo $NODE_list | cut -d ' ' -f $NODE_NUM`
    echo " >> STGIN for $n"
    ssh $n "mkdir -p $TMPDIR"

    # mkdir
    if test $NODE_NUM -eq 1 ; then
      ssh $n "mkdir -p $TMPDIR/infl_mul"
      ssh $n "mkdir -p $TMPDIR/peff_lpf"
      ssh $n "mkdir -p $TMPDIR/rsmp_lpf"
      ssh $n "mkdir -p $TMPDIR/log"
      ssh $n "mkdir -p $TMPDIR/wvec" # SK
      for MEM in mean sprd ; do
        ssh $n "mkdir -p $TMPDIR/anal/$MEM"
        ssh $n "mkdir -p $TMPDIR/anal_f/$MEM"
        ssh $n "mkdir -p $TMPDIR/gues/$MEM"
      done
      ssh $n "mkdir -p $TMPDIR/gues/kld"

      # Inflation
      if test -e $OUTPUT/infl_mul/$PYMDH.grd ; then
        echo "  > INFLATION"
        scp $OUTPUT/infl_mul/$PYMDH.grd  $n:$TMPDIR/infl_mul
      fi
      # Local PF
      if test -e $OUTPUT/peff_lpf/$PYMDH.grd ; then
        scp $OUTPUT/peff_lpf/$PYMDH.grd  $n:$TMPDIR/peff_lpf
      fi
      if test -e $OUTPUT/rsmp_lpf/$PYMDH.grd ; then
        scp $OUTPUT/rsmp_lpf/$PYMDH.grd  $n:$TMPDIR/rsmp_lpf
      fi
    fi
    # gues & inputs
    scp $TMPDIR/tmp_stgin/gues_$n.tar $n:$TMPDIR
    scp $TMPDIR/tmp_stgin/other.tar   $n:$TMPDIR
    ssh $n "cd $TMPDIR ; tar xf gues_$n.tar ; tar xf other.tar ; rm *.tar"
    CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "  STGIN($n) ended: $CTIME" >> $TIMELOG

    NODE_NUM=$((NODE_NUM + 1))
  done
  cd $TMPDIR
fi
get_ELAPS
CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "Settings ended: $CTIME, $ELAPS (s)" >> $TIMELOG

#
# Cycle run ### MAIN LOOP ========================================>>> (str) 
#
PBS_NODEFILE="./nodefile_letkf"
while test $IYYYY$IMM$IDD$IHH -le $EYYYY$EMM$EDD$EHH ; do
echo '>>>'
echo ">>> BEGIN COMPUTATION OF $IYYYY/$IMM/$IDD/$IHH"
echo '>>>'
  CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "DATE of Current run : $IYYYY/$IMM/$IDD ${IHH}:00 at $CTIME" >> $TIMELOG
  TYMDH=`date -d "$IYYYY$IMM$IDD $IHH:00 6 hour" +%Y%m%d%H`
  PYMDH=`date -d "$IYYYY$IMM$IDD $IHH:00 6 hour ago" +%Y%m%d%H`
  TY=`echo $TYMDH | cut -c1-4` ; TM=`echo $TYMDH | cut -c5-6`
  TD=`echo $TYMDH | cut -c7-8` ; TH=`echo $TYMDH | cut -c9-10`
#
# LETKF
#
if [ -e letkf.cnf ] ; then rm letkf.cnf ; fi
cat >> letkf.cnf << EOF
&letkf_param
  sigma_obs  = ${HSIG}.0d3,
  nbv        = ${MEMBER},
  pymdh      = $PYMDH,
  ymdh       = $IYYYY$IMM$IDD$IHH,
  logic_wout = .true.,
  logic_wsth = ${LOGIC_STH},
  logic_wint = ${LOGIC_INT},
  dastype    = $DASTYPE,
  resample_m = $SIRPEF,
  type_pfmtx = $RSM,
  type_relax = $RTP,
  alph_relax = $RLX,
  fgt_factor = $FGT,
  gamma_gmpf = ${LFACTOR}d0,  
/
EOF

echo "  > OBSOPE"
set_START
#./$OBSOPE > obsope.log
#mpirun -np $NODE ./$OBSOPE #> obsope.log
$MPIRUN -np $NODE ./$OBSOPE
get_ELAPS
CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "OBSOPE : $CTIME, $ELAPS (s)" >> $TIMELOG

echo "  > LETKF"
set_START
#mpirun -np $PE ./$LETKF #< /dev/null
#mpirun -np $NODE ./$LETKF #< /dev/null
time $MPIRUN -np $NODE ./$LETKF
tail -n 26 NOUT-000000
mv NOUT-000000 $TMPDIR/log/$IYYYY$IMM$IDD$IHH.log
get_ELAPS
CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "LETKF  : $CTIME, $ELAPS (s)" >> $TIMELOG

#
# ensemble forecast
#
CTIME=$(date +%Y/%m/%d-%T)
echo " >>"
echo " >> Ensemple Prediction    (SPEEDY)"
echo " >>"
echo " >>"
echo " >> ENSEMBLE PREDICTION with $NODE procs"
echo " >>"
set_START
sh 110_ensfcst_MPI.sh $IYYYY$IMM$IDD$IHH $MEMBER $NODE
get_ELAPS
CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "FCST   : $CTIME,ini $ELAPS (s)" >> $TIMELOG

echo " >>"
echo " >> COMPUTE ENSEMBLE MEAN"
echo " >>"
set_START
#if test $DISK -eq 0 ; then
#  sh 120_ensemble_mean_MPI.sh $IYYYY$IMM$IDD$IHH $MEMBER $EXP $OUTPUT $PE
#elif test $DISK -eq 1 ; then
  sh 120_ensemble_mean_MPI.sh $IYYYY$IMM$IDD$IHH $MEMBER $EXP $TMPDIR $NODE
#fi
get_ELAPS
CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "ENSmean: $CTIME, $ELAPS (s)" >> $TIMELOG
#
# Clean up
#
if test $NOSAVEENS -eq 1 ; then
  if test `echo $PYMDH | cut -c7-10` != "0100" ; then
    for n in $NODE_list ; do
      ssh $n "rm -f $TMPDIR/gues/??????/${PYMDH}.grd   $TMPDIR/gues/??????/${PYMDH}_p.grd   \
                    $TMPDIR/anal_f/??????/${PYMDH}.grd $TMPDIR/anal_f/??????/${PYMDH}_p.grd \
                    $TMPDIR/anal/??????/${PYMDH}.grd   $TMPDIR/anal/??????/${PYMDH}_w.grd" &
    done
  fi
fi # NOSAVEENS
#
# Date change ### MAIN LOOP END ###
#
IYYYY=$TY ; IMM=$TM ; IDD=$TD ; IHH=$TH
done
wait
#
# Cycle run ### MAIN LOOP ========================================>>> (end) 
#
# STAGE OUT
#
if test $DISK -eq 1 ; then
  get_ELAPS
  CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "STGOUTs: $CTIME, $ELAPS (s)" >> $TIMELOG

  echo '>>>'
  echo '>>> STAGE OUT'
  NODE_NUM=1
  for n in $NODE_list ; do
    DIRLIST_STGOUT="gues anal anal_f infl_mul log wvec peff_lpf rsmp_lpf"
    echo " >> STGOUT for $n"
    if test $NODE_NUM -eq 1 ; then
#      ssh $n "cd $TMPDIR ; tar cf stgout.tar gues anal anal_f infl_mul log wvec"
#      scp $n:$TMPDIR/stgout.tar $OUTPUT/stgout_$n.tar
      for dirname in ${DIRLIST_STGOUT} ; do
        echo " >> STGOUT for $n $dirname"
        ssh $n "cd $TMPDIR ; tar cf stgout_${dirname}.tar ${dirname}"
        scp $n:$TMPDIR/stgout_${dirname}.tar $OUTPUT/stgout_${dirname}_$n.tar
        rm stgout_${dirname}.tar
      done
#    else
#      ssh $n "cd $TMPDIR ; tar cf stgout.tar gues anal anal_f infl_mul log wvec"
#      scp $n:$TMPDIR/stgout.tar $OUTPUT/stgout_$n.tar
    fi
    NODE_NUM=$((NODE_NUM + 1))
  done

  get_ELAPS
  CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "TAROUTs: $CTIME, $ELAPS (s)" >> $TIMELOG
  cd $OUTPUT
  for n in $NODE_list ; do
#    tar xf stgout_$n.tar &
    for dirname in ${DIRLIST_STGOUT} ; do
      tar xf stgout_${dirname}_$n.tar &
    done
  done

# CLEAR MEMORY SPACE on COMPUTATINAL NODES
  echo '>>>'
  echo '>>> CLEAR MEMORY SPACE'
  for n in $NODE_list ; do
    ssh $n "rm -fr $TMPDIR" &
  done
  wait

  for n in $NODE_list ; do
    rm -f stgout_$n.tar &
  done
fi

if test $DISK -eq 0 ; then
  #rm -fr $SPEEDY/DATA/tmp_$OBS_$EXP
  echo "skip to remove tmpdir"
elif test $DISK -eq 1 ; then
  rm -fr /dev/shm/kotsuki_${OBS}_${EXP}
fi
wait


rm -fr $CDIR/nodefile/$EXP
date +"%Y/%m/%d %H:%M:%S"
CTIME=$(date +"%Y/%m/%d %H:%M:%S") ; echo "Job ended at $CTIME" >> $TIMELOG
echo "NORMAL END :: OBS & EXP :: " $OBS $EXP

####
done; #dastbl
done ; done ; done               #HSG   #OBS   #MEM
done ; done                      #SMP   #LUT
done ; done ; done ; done ; done #RTP   #RLX   #RSM   #FGT   #STH

echo "end all in jobs"
exit
