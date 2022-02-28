PROGRAM letkf
!=======================================================================
!
! [PURPOSE:] Main program of LETKF
!
! [HISTORY:]
!   01/16/2009 Takemasa Miyoshi  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_time
  USE common_speedy
  USE common_mpi_speedy
  USE common_letkf
  USE letkf_obs
  USE letkf_tools
  USE interpolate

  IMPLICIT NONE
  REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: gues2d(:,:,:)
  REAL(r_size),ALLOCATABLE :: anal3d(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: anal2d(:,:,:)
  !!REAL(r_size) :: rtimer00,rtimer
  INTEGER :: ierr
  CHARACTER(11) :: stdoutf='NOUT-000000'
  CHARACTER(4)  :: guesf='gs00'

  namelist / letkf_param / sigma_obs, nbv, pymdh, ymdh,                    &
    logic_wout, logic_wsth, logic_wint, dastype,                           &
    resample_m, type_pfmtx, type_relax, alph_relax, fgt_factor, gamma_gmpf
!-----------------------------------------------------------------------
! Initial settings
!-----------------------------------------------------------------------
!TMP!  rtimer00 = MPI_WTIME() ; rtimer01 = rtimer00
  CALL initialize_mpi
  CALL set_timer
  rtimer00 = MPI_WTIME() ; rtimer01 = rtimer00 ! SK 20200804
 
 !==> SK 20180607 for exp.
  open(1,file='letkf.cnf')
    read(1,nml=letkf_param)
  close(1)
!
  WRITE(stdoutf(9:11), '(I3.3)') myrank
  WRITE(6,'(3A,I3.3)') 'STDOUT goes to ',stdoutf,' for MYRANK ', myrank
  OPEN(6,FILE=stdoutf)
  WRITE(6,'(A,I3.3,2A)') 'MYRANK=',myrank,', STDOUTF=',stdoutf
!
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '  LOCAL ENSEMBLE TRANSFORM KALMAN FILTERING  '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '   LL      EEEEEE  TTTTTT  KK  KK  FFFFFF    '
  WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(6,'(A)') '   LL      EEEEE     TT    KKK     FFFFF     '
  WRITE(6,'(A)') '   LL      EE        TT    KK KK   FF        '
  WRITE(6,'(A)') '   LLLLLL  EEEEEE    TT    KK  KK  FF        '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '             WITHOUT LOCAL PATCH             '
  WRITE(6,'(A)') '                                             '
  WRITE(6,'(A)') '          Coded by Takemasa Miyoshi          '
  WRITE(6,'(A)') '  Based on Ott et al (2004) and Hunt (2005)  '
  WRITE(6,'(A)') '  Tested by Miyoshi and Yamane (2006)        '
  WRITE(6,'(A)') '============================================='
  WRITE(6,'(A)') '              LETKF PARAMETERS               '
  WRITE(6,'(A)') ' ------------------------------------------- '
  WRITE(6,'(A,I15)')   '   nbv             :',nbv
  WRITE(6,'(A,I15)')   '   nslots          :',nslots
  WRITE(6,'(A,I15)')   '   nbslot          :',nbslot
  WRITE(6,'(A,F15.2)') '   sigma_obs       :',sigma_obs
  WRITE(6,'(A,F15.2)') '   sigma_obsv      :',sigma_obsv
  WRITE(6,'(A,F15.2)') '   sigma_obst      :',sigma_obst
  WRITE(6,'(A,i)')     '   cymdh           :',ymdh    
  WRITE(6,'(A,i)')     '   pymdh           :',pymdh   
  WRITE(6,'(A,f)')     '   sigma_obs       :',sigma_obs
  WRITE(6,'(A,L)')     '   logic_wout      :',logic_wout
  WRITE(6,'(A,L)')     '   logic_wsth      :',logic_wsth
  WRITE(6,'(A,L)')     '   logic_wint      :',logic_wint
  WRITE(6,'(A,I,A)')   '   type_relax      :',type_relax, "  0:no, 1:RTPP, 2:RTPS"
  WRITE(6,'(A,F)')     '   alph_relax      :',alph_relax

  WRITE(6,'(A,I)')     '   LPF,resample_m  :',resample_m
  WRITE(6,'(A,F)')     '   LPF,fgt_factor  :',fgt_factor
  WRITE(6,'(A,F)')     '   LPF,gamma_gmpf  :',gamma_gmpf    
  WRITE(6,'(A,I,A)')   '   LPF,type_pfmtx  :',type_pfmtx, "  1: Potthast et al. (2019;MWR)  2: Kondo et al. (2020)"
  WRITE(6,'(A,I)')     '   LPF,dastype     :',dastype
  IF( dastype==0 ) WRITE(6,'(A)') '   dastype=0 ::   mean update :: LETKF,   ptb update :: LETKF'
  IF( dastype==1 ) WRITE(6,'(A)') '   dastype=1 ::   mean update :: LETKF,   ptb update :: DEBUG, LPF only for mean update'
  IF( dastype==2 ) WRITE(6,'(A)') '   dastype=2 ::   mean update :: LPF  ,   ptb update :: LAPF  resampling'
  IF( dastype==3 ) WRITE(6,'(A)') '   dastype=3 ::   mean update :: LPF  ,   ptb update :: LPFGM resampling'

  WRITE(6,'(A)') '============================================='
  CALL set_common_speedy
  CALL set_common_mpi_speedy
  ALLOCATE(gues3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(gues2d(nij1,nbv,nv2d))
  ALLOCATE(anal3d(nij1,nlev,nbv,nv3d))
  ALLOCATE(anal2d(nij1,nbv,nv2d))
!
  rtimer = MPI_WTIME()
!  WRITE(6,'(A,2F10.2)') '### TIMER(INITIALIZE):',rtimer-rtimer01,rtimer-rtimer00
  rtimerl(1) = rtimer-rtimer01
  rtimer01 = rtimer
!-----------------------------------------------------------------------
! Observations
!-----------------------------------------------------------------------
  !
  ! CONVENTIONAL OBS
  !
  CALL set_letkf_obs
!
  rtimer = MPI_WTIME()
!  WRITE(6,'(A,2F10.2)') '### TIMER(READ_OBS):',rtimer-rtimer01,rtimer-rtimer00
  rtimerl(2) = rtimer-rtimer01
  rtimer01=rtimer
!-----------------------------------------------------------------------
! First guess ensemble
!-----------------------------------------------------------------------
  !
  ! READ GUES
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  WRITE(guesf(3:4),'(I2.2)') nbslot
  CALL read_ens_mpi(guesf,nbv,gues3d,gues2d)
  !
  ! WRITE ENS MEAN and SPRD
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL write_ensmspr_mpi('gues',nbv,gues3d,gues2d)
!
  rtimer = MPI_WTIME()
!  WRITE(6,'(A,2F10.2)') '### TIMER(READ_GUES):',rtimer-rtimer01,rtimer-rtimer00
  rtimerl(3) = rtimer-rtimer01
  rtimer01=rtimer
!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
  !
  ! LETKF
  !

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
! --------------->
  WRITE(6,*) 'letkf begin'
  CALL das_letkf(gues3d,gues2d,anal3d,anal2d)
! --------------->
  WRITE(6,*) 'letkf done'
  rtimer = MPI_WTIME()
!  WRITE(6,'(A,2F10.2)') '### TIMER(DAS_LETKF):',rtimer-rtimer01,rtimer-rtimer00
  rtimerl(4) = rtimer-rtimer01
  rtimer01=rtimer
!-----------------------------------------------------------------------
! Analysis ensemble
!-----------------------------------------------------------------------
  !
  ! WRITE ANAL
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL write_ens_mpi('anal',nbv,anal3d,anal2d)
  !
  ! WRITE ENS MEAN and SPRD
  !
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL write_ensmspr_mpi('anal',nbv,anal3d,anal2d)
!
  rtimer = MPI_WTIME()
!  WRITE(6,'(A,2F10.2)') '### TIMER(WRITE_ANAL):',rtimer-rtimer01,rtimer-rtimer00
  rtimerl(5) = rtimer-rtimer01
  rtimer01=rtimer
!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
  !SK follows KK!CALL monit_mean('gues')
  !SK follows KK!CALL monit_mean('anal')
!
  rtimer = MPI_WTIME()
! WRITE(6,'(A,2F10.2)') '### TIMER(MONIT_MEAN):',rtimer-rtimer01,rtimer-rtimer00
  rtimerl(6) = rtimer-rtimer01
  rtimer01=rtimer

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_Gather(rtimerl,slot,MPI_REAL8,rtimer_all,slot,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  !KK!CALL MPI_Gather(ltimerl,slot+1,MPI_REAL8,ltimer_all,slot+1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_00)  :', ltimerl(0),  ltimerl(0)/ sum(rtimerl(:))*100, ' %', ltimerl(0)/ sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_01)  :', ltimerl(1),  ltimerl(1)/ sum(rtimerl(:))*100, ' %', ltimerl(1)/ sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_02)  :', ltimerl(2),  ltimerl(2)/ sum(rtimerl(:))*100, ' %', ltimerl(2)/ sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_03)  :', ltimerl(3),  ltimerl(3)/ sum(rtimerl(:))*100, ' %', ltimerl(3)/ sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_04)  :', ltimerl(4),  ltimerl(4)/ sum(rtimerl(:))*100, ' %', ltimerl(4)/ sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_05)  :', ltimerl(5),  ltimerl(5)/ sum(rtimerl(:))*100, ' %', ltimerl(5)/ sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_06)  :', ltimerl(6),  ltimerl(6)/ sum(rtimerl(:))*100, ' %', ltimerl(6)/ sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_07)  :', ltimerl(7),  ltimerl(7)/ sum(rtimerl(:))*100, ' %', ltimerl(7)/ sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_08)  :', ltimerl(8),  ltimerl(8)/ sum(rtimerl(:))*100, ' %', ltimerl(8)/ sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_09)  :', ltimerl(9),  ltimerl(9)/ sum(rtimerl(:))*100, ' %', ltimerl(9)/ sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_10)  :', ltimerl(10), ltimerl(10)/sum(rtimerl(:))*100, ' %', ltimerl(10)/sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_11)  :', ltimerl(11), ltimerl(11)/sum(rtimerl(:))*100, ' %', ltimerl(11)/sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_12)  :', ltimerl(12), ltimerl(12)/sum(rtimerl(:))*100, ' %', ltimerl(12)/sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_13)  :', ltimerl(13), ltimerl(13)/sum(rtimerl(:))*100, ' %', ltimerl(13)/sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_14)  :', ltimerl(14), ltimerl(14)/sum(rtimerl(:))*100, ' %', ltimerl(14)/sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_15)  :', ltimerl(15), ltimerl(15)/sum(rtimerl(:))*100, ' %', ltimerl(15)/sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkf_16)  :', ltimerl(16), ltimerl(16)/sum(rtimerl(:))*100, ' %', ltimerl(16)/sum(ltimerl(:))*100, ' %'
!KK!    WRITE(6,'(A,2F10.2,A,F10.2,A)') '### TIMER(letkfsum)  :', sum(ltimerl(:)), sum(ltimerl(:))/sum(rtimerl(:))*100, ' %', sum(ltimerl(:))/sum(ltimerl(:))*100, ' %'

    WRITE(6,*) ''
    WRITE(6,'(A,2F10.2,A)') '### TIMER(INITIALIZE):', rtimerl(1), rtimerl(1)/sum(rtimerl(:))*100, ' %'
    WRITE(6,'(A,2F10.2,A)') '### TIMER(READ_OBS)  :', rtimerl(2), rtimerl(2)/sum(rtimerl(:))*100, ' %'
    WRITE(6,'(A,2F10.2,A)') '### TIMER(READ_GUES) :', rtimerl(3), rtimerl(3)/sum(rtimerl(:))*100, ' %'
    WRITE(6,'(A,2F10.2,A)') '### TIMER(DAS_LETKF) :', rtimerl(4), rtimerl(4)/sum(rtimerl(:))*100, ' %'
    WRITE(6,'(A,2F10.2,A)') '### TIMER(WRITE_ANAL):', rtimerl(5), rtimerl(5)/sum(rtimerl(:))*100, ' %'
    WRITE(6,'(A,2F10.2,A)') '### TIMER(MONIT_MEAN):', rtimerl(6), rtimerl(6)/sum(rtimerl(:))*100, ' %'
    WRITE(6,'(A,2F10.2,A)') '### TIMER(TOTAL)     :', sum(rtimerl(:)), sum(rtimerl(:))/sum(rtimerl(:))*100, ' %'

!-----------------------------------------------------------------------
! Finalize
!-----------------------------------------------------------------------
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi
  STOP
END PROGRAM letkf
