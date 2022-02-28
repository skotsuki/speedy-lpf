#!/bin/sh
#==================================================#
export ONLYCAT=1 # 0: no  1: only catting mode
#==================================================#

set -e
CWD=`pwd`
cd ../ ; SPEEDY=`pwd`
cd $CWD
if [ ! -d ./grddata ] ; then mkdir -p ./grddata ; fi

CATDIR=$CWD/gnu_cat ; if [ ! -d $CATDIR ] ; then mkdir -p $CATDIR ; fi
OUTDIR=$CWD/gnu_000 ; if [ ! -d $OUTDIR ] ; then mkdir -p $OUTDIR ; fi

# ========>>> preparation           <<<============#
FOPT='-convert big_endian -assume byterecl -O3'
###FOPT='-convert big_endian -assume byterecl -fopenmp'
if [ -e prg_gues-rmse.exe ] ; then rm prg_gues-rmse.exe ; fi
ifort $FOPT prg_gues-rmse.f90 -o prg_gues-rmse.exe

# ========>>> parameter for spin-up <<<============#
ADATE=1982010100 # SET FOR SPIN-UP
SPNDIR=$SPEEDY/DATA/reg2/BASELETKF_M000040L0500IADP
CTLDIR=$SPEEDY/DATA/reg2/BASELETKF_M000040L0500IADP

# ========>>> parameter for exps    <<<============#
NRUNME=28         # DALYED MEAN
SDATE=1982020100 # SET FOR EXPERIMENT (STR)
EDATE=1982053118 # SET FOR EXPERIMENT (END)
#EDATE=1982033118 # SET FOR EXPERIMENT (END)

PDATE=1982030100 # SET FOR AVE RMSE   (STR)
QDATE=$EDATE     # SET FOP AVE RMSE   (END)

# ========>>> parameter for basics  <<<============#
MEMTBL="40" ; OBSTBL="raob" ; HSGTBL="800"
MEMTBL="40" ; OBSTBL="reg2" ; HSGTBL="500"

# ========>>> parameter for LPFGM   <<<============#
SYS="SYS20200418_" # parameter sweeep
SYS="SYS20200425_" # resampling matrix revised (KK)
SYS="SYS20200426_" # ramdum number fixed
SYS="SYS20200427_" # weight-succession revised, experiments

FACTBL="2.5"
WSMTBL="woWSM"
DASTBL="3LPFGM 2LAPF_"
RTPTBL="RTPS"

#====> LPFGM TABLE
DASTBL="3LPFGM"
FGTTBL="1.00 0.00 0.80"
INFTBL="infKK"
NEFTBL="2 5 10 40"
ALPTBL="0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00"

#====> LAPF  TABLE
#DASTBL="2LAPF_"
#HSGTBL="500 400 300 200"
#FGTTBL="0.00"
#INFTBL="infKK"
#NEFTBL="40"
#ALPTBL="0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.15 1.20"

##====> debug TABLE
#DASTBL="1LPFck"
#HSGTBL="500 400 300"
#FGTTBL="0.00"
#INFTBL="infKK"
#NEFTBL="40"
#ALPTBL="0.00"

############################################################################ (1) BASIC LOOPS
for MEMBER in $MEMTBL ; do
for OBS    in $OBSTBL ; do
for HSG    in $HSGTBL ; do
  M=$( printf %06i ${MEMBER} )
  H=$( printf %04i ${HSG} )
  BASIC=M${M}L${H}IADP
  BASIC_LETKF=M${M}L0500IADP
############################################################################ (1) BASIC (calc. rmses always for CTRL)
cd $CWD
  OUTTIME=${OUTDIR}/time_${OBS}_LETKF_${BASIC_LETKF}.txt
  OUTTAVE=${OUTDIR}/tave_${OBS}_LETKF_${BASIC_LETKF}_SMP${PDATE}-${QDATE}.txt
  N=$( printf %03i ${MEMBER} )
if [ -e speedy.cnf ] ; then rm speedy.cnf ; fi
cat >> speedy.cnf << EOF
&speedy_param
  adate         = $ADATE,
  sdate         = $SDATE,
  edate         = $EDATE,
  pdate         = $PDATE,
  qdate         = $QDATE,
  nbv           = $MEMBER,
  nrunme        = $NRUNME, 
  natdir        = '$SPEEDY/DATA/nature/'
  spndir        = '${SPNDIR}',
  expdir        = '${CTLDIR}',
  outtime       = '${OUTTIME}',
  outtave       = '${OUTTAVE}',
  !expname       = '${EXPNAME}',
  expname       = 'LETKF_${BASIC_LETKF}',
  obsname       = '${OBS}'
  cneff         = '$N',
  !crtps         = '$ALP',
  crtps         = '0.00',
/
EOF
if [ $ONLYCAT -ne 1 ] ; then  ./prg_gues-rmse.exe ; fi
############################################################################ (2) EXPERIMENTAL LOOPS (DAS)
for DAS in $DASTBL ; do
for RTP in $RTPTBL ; do
for WSM in $WSMTBL ; do
############################################################################ (3) EXPERIMENTAL LOOPS (LPF)
for FAC in $FACTBL ; do
for FGT in $FGTTBL ; do
for INF in $INFTBL ; do
for NEF in $NEFTBL ; do
  N=$( printf %03i ${NEF} )
  CATNAME=${DAS}_${BASIC}_${RTP}x.xx_${WSM}_${INF}_fac${FAC}_fgt${FGT}_rsmp${N}
  if [ $DAS = "1LPFck" ]; then
  CATNAME=${DAS}_${BASIC}_${RTP}x.xx_${WSM}                                     ; fi
    
  CATFILE=${CATDIR}/tave_${OBS}_${CATNAME}_SMP${PDATE}-${QDATE}.txt
  if [ -e $CATFILE ] ; then rm $CATFILE ; fi
  touch $CATFILE
for ALP in $ALPTBL ; do
  EXPNAME=${DAS}_${BASIC}_${RTP}${ALP}_${WSM}_${INF}_fac${FAC}_fgt${FGT}_rsmp${N}
  if [ $DAS = "1LPFck" ]; then
  EXPNAME=${DAS}_${BASIC}_${RTP}${ALP}_${WSM}                                     ; fi

  EXPDIR=$SPEEDY/DATA/$OBS/$SYS$EXPNAME

  OUTTIME=${OUTDIR}/time_${OBS}_${EXPNAME}.txt
  OUTTAVE=${OUTDIR}/tave_${OBS}_${EXPNAME}_SMP${PDATE}-${QDATE}.txt

cd $CWD
if [ -e speedy.cnf ] ; then rm speedy.cnf ; fi
cat >> speedy.cnf << EOF
&speedy_param
  adate         = $ADATE,
  sdate         = $SDATE,
  edate         = $EDATE,
  pdate         = $PDATE,
  qdate         = $QDATE,
  nbv           = $MEMBER,
  nrunme        = $NRUNME, 
  natdir        = '$SPEEDY/DATA/nature/'
  spndir        = '${SPNDIR}',
  expdir        = '${EXPDIR}',
  outtime       = '${OUTTIME}',
  outtave       = '${OUTTAVE}',
  expname       = '${EXPNAME}',
  obsname       = '${OBS}'
  cneff         = '$N',
  crtps         = '$ALP',
/
EOF
if [ $ONLYCAT -ne 1 ] ; then  ./prg_gues-rmse.exe                             ; fi
  cat ${OUTTAVE} >> $CATFILE
if [ $ONLYCAT -eq 1 ] ; then  echo "  catting " ${OUTTAVE} " into " $CATFILE  ; fi

done                        # ALP
done ; done ; done ; done   # FAC, FGT, INF, NEF
done ; done ; done          # DAS, RTP, WSM
done ; done ; done          # MEM, OBS, HSG
