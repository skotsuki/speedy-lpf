MODULE letkf_tools
!=======================================================================
!
! [PURPOSE:] Module for LETKF with SPEEDY
!
! [HISTORY:]
!   01/26/2009 Takemasa Miyoshi  created
!
!=======================================================================
  USE common
  USE common_mpi
  USE common_time
  USE common_speedy
  USE common_mpi_speedy
  USE common_letkf
  USE common_lpf
  USE letkf_obs
  USE interpolate

  IMPLICIT NONE

  PRIVATE
  PUBLIC ::  das_letkf

  INTEGER,SAVE :: nobstotal

  REAL(r_size),PARAMETER :: cov_infl_mul = -1.01d0 !multiplicative inflation
! > 0: globally constant covariance inflation
! < 0: 3D inflation values input from a GPV file "infl_mul.grd"
  REAL(r_size),PARAMETER :: sp_infl_add = 0.d0 !additive inflation
!TVS  LOGICAL,PARAMETER :: msw_vbc = .FALSE.
  REAL(r_size),PARAMETER :: var_local(nv3d+nv2d,nid_obs) = RESHAPE( &
!           U      V      T      Q     PS   RAIN
   & (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  & ! U
   &    1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  & ! V
   &    1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  & ! T
   &    1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  & ! Q
   &    1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  & ! RH
   &    1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,  & ! PS
   &    1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)& ! RAIN
   & ,(/nv3d+nv2d,nid_obs/))
  INTEGER,SAVE :: var_local_n2n(nv3d+nv2d)

CONTAINS
!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
SUBROUTINE das_letkf(gues3d,gues2d,anal3d,anal2d)
  IMPLICIT NONE
  CHARACTER(12) :: inflfile='infl_mul.grd'
  CHARACTER(29) :: wvecfile='./wvec/yyyymmddhh_MXXXXXX.grd'             !-->A.P 6/6/2018 weight vec grd file
  CHARACTER(25) :: inflfile_i='./infl_mul/yyyymmddhh.grd'
  CHARACTER(25) :: inflfile_o='./infl_mul/yyyymmddhh.grd'
  REAL(r_size),INTENT(INOUT) :: gues3d(nij1,nlev,nbv,nv3d) ! background ensemble
  REAL(r_size),INTENT(INOUT) :: gues2d(nij1,nbv,nv2d)      !  output: destroyed
  REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,nbv,nv3d)   ! analysis ensemble
  REAL(r_size),INTENT(OUT) :: anal2d(nij1,nbv,nv2d)
  REAL(r_size),ALLOCATABLE :: mean3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: mean2d(:,:)
  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:)
  REAL(r_size),ALLOCATABLE :: rloc(:)
  REAL(r_size),ALLOCATABLE :: dep(:)
  REAL(r_size),ALLOCATABLE :: work3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: work2d(:,:)
  REAL(r_sngl),ALLOCATABLE :: work3dg(:,:,:,:)
  REAL(r_sngl),ALLOCATABLE :: work2dg(:,:,:)
  REAL(r_size),ALLOCATABLE :: logpfm(:,:)
  REAL(r_size) :: parm
  REAL(r_size) :: trans(nbv,nbv,nv3d+nv2d)
  LOGICAL :: ex
  INTEGER :: ij,ilev,n,m,i,j,k,nobsl,ierr


  !==> LPF ; SK 20190420
  LOGICAL      :: logic_lpfgm = .false.
  INTEGER      :: iseed
  REAL(r_size) :: srand(nbv)
  REAL(r_size) :: gr_common(nbv,nbv)            , Ident(nbv,nbv)
!STACK!  REAL(r_size) :: nobs3d(nij1,nlev,nv3d)        , nobs2d(nij1,nv2d)
!STACK!  REAL(r_size) :: wvec3d(nij1,nlev,nv3d,nbv)    , wvec2d(nij1,nv2d,nbv)
!STACK!  REAL(r_size) :: wmat3d(nij1,nlev,nv3d,nbv,nbv), wmat2d(nij1,nv2d,nbv,nbv)
!STACK!  REAL(r_size) :: pvec3d(nij1,nlev,nv3d,nbv)    , pvec2d(nij1,nv2d,nbv)    
!STACK!  REAL(r_size) :: pmat3d(nij1,nlev,nv3d,nbv,nbv), pmat2d(nij1,nv2d,nbv,nbv)
  REAL(r_size), ALLOCATABLE :: nobs3d(:,:,:)        , nobs2d(:,:)
  REAL(r_size), ALLOCATABLE :: wvec3d(:,:,:,:)    , wvec2d(:,:,:)
  REAL(r_size), ALLOCATABLE :: wmat3d(:,:,:,:,:), wmat2d(:,:,:,:)
  REAL(r_size), ALLOCATABLE :: pvec3d(:,:,:,:)    , pvec2d(:,:,:)    
  REAL(r_size), ALLOCATABLE :: pmat3d(:,:,:,:,:), pmat2d(:,:,:,:)

  REAL(r_size) :: asis3d(nij1,nlev,nv3d,nbv)    , asis2d(nij1,nv2d,nbv)  
  REAL(r_size) :: peff3d(nij1,nlev,nv3d)        , peff2d(nij1,nv2d)
  REAL(r_size) :: rsmp3d(nij1,nlev,nv3d)        , rsmp2d(nij1,nv2d)
  CHARACTER(25) :: pefffile_o='./peff_lpf/yyyymmddhh.grd'
  CHARACTER(25) :: rsmpfile_o='./rsmp_lpf/yyyymmddhh.grd'
  CHARACTER(33) :: asisfile_o='./asis_lpf/yyyymmddhh_MXXXXXX.grd'

  !==> Weight Smoother
  REAL(r_size) :: sigma_g, distg, dist_zerog, sfnc, swgh
  REAL(r_size) :: wgh_nij2map(nij1,nlon,nlat)

  !==> Weight Interpolation
  INTEGER              :: iprocs, iexe, nexe, ilon, ilat
  INTEGER, ALLOCATABLE :: proc_m(:,:)
      
  REAL(r_size) :: trans3(nbv,nbv)
  REAL(r_size) :: trans2(nbv,nbv)
  REAL(r_size) :: msk(nlon,nlat),   msk_me(nij1) 
  REAL(r_size) :: wix(nlon,nlat,4), wix_me(nij1,4)
  REAL(r_size) :: wiy(nlon,nlat,4), wiy_me(nij1,4)
  REAL(r_size) :: fac(nlon,nlat,4), fac_me(nij1,4)
  REAL(r_size) :: tmpave, tmpinf, tmpptb(nbv), gusspr, anlspr

!
  ltimer00 = MPI_WTIME() ; ltimer01 = ltimer00
  ptimer00 = MPI_WTIME() ; ptimer01 = ptimer00
  ptimer   = ptimer00
!
  WRITE(6,'(A)') 'Hello from das_letkf'
  nobstotal = nobs !+ ntvs
  WRITE(6,'(A,I8)') 'Target observation numbers : NOBS=',nobs!,', NTVS=',ntvs
  !
  ! identity matrix
  !
  Ident(:,:) = 0.0d0
  DO k=1,nbv ; Ident(k,k) = 1.0d0 ; END DO
  !SK,NO-WI;20200417!!
  !SK,NO-WI;20200417!! gauss weight preparation for filtering
  !SK,NO-WI;20200417!!
  !SK,NO-WI;20200417!IF( logic_wsth )THEN
  !SK,NO-WI;20200417!  sigma_g            = sigma_obs * 1.0d0
  !SK,NO-WI;20200417!  dist_zerog         = sigma_g * SQRT(10.0d0/3.0d0) * 2.0d0
  !SK,NO-WI;20200417!  wgh_nij2map(:,:,:) = 0.0d0
  !SK,NO-WI;20200417!  DO ij=1,nij1 ; DO j=1,nlat ; DO i=1,nlon
  !SK,NO-WI;20200417!    CALL com_distll_1(lon1(ij),lat1(ij),lon(i),lat(j),distg)
  !SK,NO-WI;20200417!    IF( distg < dist_zerog ) &
  !SK,NO-WI;20200417!      wgh_nij2map(ij,i,j) = EXP(-0.5d0 * ((distg/sigma_g)**2) )
  !SK,NO-WI;20200417!  END DO       ; END DO      ; END DO
  !SK,NO-WI;20200417!END IF
  !
  ! preparation for member parallerization
  !
  !===> setting for parallel computing
  IF( MOD(nbv,nprocs) /= 0 ) THEN
    print *, "error, nbv should be devided by the nprocs"
    STOP ; ENDIF

  nexe = nbv / nprocs
  allocate ( proc_m(0:nprocs-1,nexe) ) ; proc_m(:,:) = -999

  iprocs = -1 ; iexe  = 1
  DO m=1,nbv
    iprocs = iprocs + 1
    IF( iprocs == nprocs ) THEN
      iprocs = 0
      iexe   = iexe + 1
    ENDIF
    proc_m(iprocs,iexe)    = m
  END DO
  !!DO iexe=1,nexe
  !!  PRINT '(3i)', iexe, myrank,proc_m(myrank,iexe)
  !!END DO

  !
  ! In case of no obs
  !
  IF(nobstotal == 0) THEN
    WRITE(6,'(A)') 'No observation assimilated'
    anal3d = gues3d
    anal2d = gues2d
    RETURN
  END IF
  !
  ! Variable localization
  !
  var_local_n2n(1) = 1
  DO n=2,nv3d+nv2d
    DO i=1,n
      var_local_n2n(n) = i
      IF(MAXVAL(ABS(var_local(i,:)-var_local(n,:))) < TINY(var_local)) EXIT
    END DO
  END DO
!print *,var_local_n2n
  !
  ! FCST PERTURBATIONS
  !
  ALLOCATE(mean3d(nij1,nlev,nv3d))
  ALLOCATE(mean2d(nij1,nv2d))
  CALL ensmean_grd(nbv,nij1,gues3d,gues2d,mean3d,mean2d)
  DO n=1,nv3d
    DO m=1,nbv
      DO k=1,nlev
        DO i=1,nij1
          gues3d(i,k,m,n) = gues3d(i,k,m,n) - mean3d(i,k,n)
        END DO
      END DO
    END DO
  END DO
  DO n=1,nv2d
    DO m=1,nbv
      DO i=1,nij1
        gues2d(i,m,n) = gues2d(i,m,n) - mean2d(i,n)
      END DO
    END DO
  END DO
  !
  ! random noises for resampling (Roland et al., 2018; MWR)
  !
  IF( dastype >= 1 .and. type_pfmtx==0 )THEN
    gr_common(:,:) = 0.0d0
    IF( myrank == 0 )THEN
      DO k=1,nbv
        CALL com_randn(nbv,gr_common(:,k))
      END DO
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST (gr_common,nbv*nbv,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr) !! share random numbers for resampling
  END IF
  !
  ! likelihood succession for local particle filter
  !
  asis3d = 1.0d0/dble(nbv)  ;  asis2d = 1.0d0/dble(nbv)
  IF( fgt_factor<1.0d0 )THEN
    !==> each member for SIS (need to save all member)
    WRITE(wvecfile,'(a7,i10.10,a2,i6.6,a4)') './asis_lpf/',pymdh,'_M',1,'.grd'
    INQUIRE(FILE=wvecfile,EXIST=ex)
    IF( ex )THEN
      ALLOCATE( work3dg(nlon,nlat,nlev,nv3d), work3d(nij1,nlev,nv3d) )
      ALLOCATE( work2dg(nlon,nlat,nv2d)     , work2d(nij1,nv2d)      )
      DO m=1,nbv
        IF(myrank == 0) THEN
          WRITE(wvecfile,'(a7,i10.10,a2,i6.6,a4)') './asis_lpf/',pymdh,'_M',m,'.grd'
          CALL read_grd4_1atm(wvecfile,work3dg,work2dg)
          IF( m==1 .or. m==nbv ) &
          WRITE(6,'(A,I3.3,3A,4f8.3)') 'MYRANK ',myrank,' is reading.. ',wvecfile, "min-max(3d,2d) :: ", &
            minval( work3dg(:,:,:,:) ), maxval( work3dg(:,:,:,:) ), &
            minval( work2dg(:,:,:)   ), maxval( work2dg(:,:,:)   )
        END IF
        CALL scatter_grd_mpi(0,work3dg,work2dg,work3d,work2d)
        asis3d(:,:,:,m) = work3d(:,:,:)
        asis2d(:,:,m)   = work2d(:,:)
      END DO
      DEALLOCATE( work3dg, work3d )
      DEALLOCATE( work2dg, work2d )

      asis3d(:,:,:,:) = asis3d(:,:,:,:)*( 1.0d0-fgt_factor ) +  fgt_factor/dble(nbv) ! forgetting-factor for likelihood
      asis2d(:,:,:)   = asis2d(:,:,:)  *( 1.0d0-fgt_factor ) +  fgt_factor/dble(nbv) ! forgetting-factor for likelihood
    ELSE
      WRITE(6,'(2A)') '!!WARNING: no such file exist: ',wvecfile
      asis3d = 1.0d0/dble(nbv)  ;  asis2d = 1.0d0/dble(nbv)
    END IF   
  ENDIF
  !
  ! multiplicative inflation
  !
  IF(cov_infl_mul > 0.0d0) THEN ! fixed multiplicative inflation parameter
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    work3d = cov_infl_mul
    work2d = cov_infl_mul
    work3d(:,nlev,:) = 1.01d0
  END IF
  IF(cov_infl_mul <= 0.0d0) THEN ! 3D parameter values are read-in
    ALLOCATE( work3dg(nlon,nlat,nlev,nv3d) )
    ALLOCATE( work2dg(nlon,nlat,nv2d) )
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    write(inflfile_i(12:21), '(I10.10)') pymdh
    INQUIRE(FILE=inflfile_i,EXIST=ex)
    IF(ex) THEN
      IF(myrank == 0) THEN
        WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',inflfile_i
        CALL read_grd4(inflfile_i,work3dg,work2dg)
        PRINT *, "  CHECK INFLATIONS :: ", minval( work3dg(:,:,:,:) ), maxval( work3dg(:,:,:,:) )
      END IF
      CALL scatter_grd_mpi(0,work3dg,work2dg,work3d,work2d)
    ELSE
      WRITE(6,'(2A)') '!!WARNING: no such file exist: ',inflfile_i
      work3d = -1.0d0 * cov_infl_mul
      work2d = -1.0d0 * cov_infl_mul
    END IF
  END IF

  !==> SK 20181015 WI based on LUT
  IF( logic_wint ) THEN
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. BILINLUT.bin for WIP'
    OPEN(1,file='BILINLUT.bin',form="unformatted",access="direct",recl=nlon*nlat*8,action="read")
      READ(1,rec=1) msk(:,:)  
      READ(1,rec=2) wix(:,:,1)  ;  READ(1,rec=6) wiy(:,:,1)  ;  READ(1,rec=10) fac(:,:,1)
      READ(1,rec=3) wix(:,:,2)  ;  READ(1,rec=7) wiy(:,:,2)  ;  READ(1,rec=11) fac(:,:,2)
      READ(1,rec=4) wix(:,:,3)  ;  READ(1,rec=8) wiy(:,:,3)  ;  READ(1,rec=12) fac(:,:,3)
      READ(1,rec=5) wix(:,:,4)  ;  READ(1,rec=9) wiy(:,:,4)  ;  READ(1,rec=13) fac(:,:,4)
    CLOSE(1)

    msk_me(:)=UNDEF ; wix_me(:,:)=UNDEF ; wiy_me(:,:)=UNDEF ; fac_me(:,:)=UNDEF
    DO m=1,nprocs
      DO i=1,nij1node(m)
        j    = m-1 + nprocs * (i-1)
        ilon = MOD(j,nlon) + 1
        ilat = (j-ilon+1) / nlon + 1
        IF( myrank==m-1 ) THEN
          msk_me(i)     = msk(ilon,ilat)
          wix_me(i,1:4) = wix(ilon,ilat,1:4)
          wiy_me(i,1:4) = wiy(ilon,ilat,1:4)
          fac_me(i,1:4) = fac(ilon,ilat,1:4)
        END IF
      END DO
    END DO
  END IF

  !
  ! p_full for background ensemble mean
  !
  ALLOCATE(logpfm(nij1,nlev))
  CALL calc_pfull(nij1,1,mean2d(:,iv2d_ps),logpfm)
  logpfm = DLOG(logpfm)
!
  ptimer = MPI_WTIME()
  WRITE(6,'(A,2F10.2)') '### TIMER in DAS_LETKF (INITIALIZE):',ptimer-ptimer01,ptimer-ptimer00
  ptimer01 = ptimer
  !
  ! MAIN ASSIMILATION LOOP
  !
  IF( dastype==3 ) logic_lpfgm = .true.
  rsmp3d(1:nij1,1:nlev,1:nv3d) = 1.0d0 ! resampling (default) 
  rsmp2d(1:nij1,       1:nv2d) = 1.0d0 ! resampling (default)

  !
  ! To fix random number for LPF
  !
  call mk_iseed(ymdh, 1,1,1,iseed)
  call com_rand_seed(nbv,iseed,srand)
  !!!print *, "inp",myrank, iseed
  
  !
  ! mtx allocation
  !
  ALLOCATE( nobs3d(nij1,nlev,nv3d)        , nobs2d(nij1,nv2d)         )
  ALLOCATE( wvec3d(nij1,nlev,nv3d,nbv)    , wvec2d(nij1,nv2d,nbv)     )
  ALLOCATE( wmat3d(nij1,nlev,nv3d,nbv,nbv), wmat2d(nij1,nv2d,nbv,nbv) )
  ALLOCATE( pvec3d(nij1,nlev,nv3d,nbv)    , pvec2d(nij1,nv2d,nbv)     )
  ALLOCATE( pmat3d(nij1,nlev,nv3d,nbv,nbv), pmat2d(nij1,nv2d,nbv,nbv) )

  ALLOCATE( hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
  DO ilev=1,nlev
    !!WRITE(6,'(A,I3)') 'ilev = ',ilev
    DO ij=1,nij1
      DO n=1,nv3d


        IF(var_local_n2n(n) < n) THEN
          trans(:,:,n)                 = trans(:,:,var_local_n2n(n))
          work3d(ij,ilev,n)            = work3d(ij,ilev,var_local_n2n(n))
          wvec3d(ij,ilev,n,1:nbv)      = wvec3d(ij,ilev,var_local_n2n(n),1:nbv)          
          wmat3d(ij,ilev,n,1:nbv,1:nbv)= wmat3d(ij,ilev,var_local_n2n(n),1:nbv,1:nbv)    
          nobs3d(ij,ilev,n)            = nobs3d(ij,ilev,var_local_n2n(n))               
          
          pvec3d(ij,ilev,n,1:nbv)      = pvec3d(ij,ilev,var_local_n2n(n),1:nbv)
          pmat3d(ij,ilev,n,1:nbv,1:nbv)= pmat3d(ij,ilev,var_local_n2n(n),1:nbv,1:nbv)
          peff3d(ij,ilev,n)            = peff3d(ij,ilev,var_local_n2n(n))
        ELSE
          CALL obs_local(ij,ilev,n,hdxf,rdiag,rloc,dep,nobsl,logpfm)
          
          !work ; LETKG (dastype=0) :: adaptive covariance inflation of letkf_core
          !work ; DEBUG (dastype=1) :: adaptive covariance inflation of letkf_core
          !work ; LAPF  (dastype=2) :: adaptive resampling amplitude of lpf_core   (used if type_pfmtx==0 i.e., LAPF's resamplig mtx)
          !work ; LPFGM (dastype=3) :: adaptive resampling amplitude of lpf_core   (used if type_pfmtx==0 i.e., LAPF's resamplig mtx)
          parm = work3d(ij,ilev,n)
          IF( dastype==3 ) parm = GAMMA_GMPF     !! LPFGM ; Walter and Potthast (2019)

          nobs3d(ij,ilev,n) = REAL( nobsl, r_size )
          CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,n),logic_lpfgm, &
                          wvec3d(ij,ilev,n,1:nbv),wmat3d(ij,ilev,n,1:nbv,1:nbv)              ) !<--- 6/6/2018 A.P. Added savestate for weight LETKF mat/vec
          IF( dastype==0 .or. dastype==1 ) work3d(ij,ilev,n) = parm !! LETKF, UPDATE 

          IF( dastype==2 .or. dastype==3 ) parm = work3d(ij,ilev,n) ! LAPF or LPFGM
          CALL lpf_core  (nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,gr_common,type_pfmtx, &
                          asis3d(ij,ilev,n,1:nbv),peff3d(ij,ilev,n),pvec3d(ij,ilev,n,1:nbv),pmat3d(ij,ilev,n,1:nbv,1:nbv))          
          IF( dastype==2 .or. dastype==3 ) work3d(ij,ilev,n) = parm ! LAPF or LPFGM
        END IF        
      END DO
      IF(ilev == 1) THEN !update 2d variable at ilev=1
        DO n=1,nv2d
          IF(var_local_n2n(nv3d+n) <= nv3d) THEN
            trans(:,:,nv3d+n)        = trans(:,:,var_local_n2n(nv3d+n))
            work2d(ij,n)             = work3d(ij,ilev,var_local_n2n(nv3d+n))
            wvec2d(ij,n,1:nbv)       = wvec3d(ij,ilev,var_local_n2n(nv3d+n),1:nbv)             
            wmat2d(ij,n,1:nbv,1:nbv) = wmat3d(ij,ilev,var_local_n2n(nv3d+n),1:nbv,1:nbv)       
            nobs2d(ij,n)             = nobs3d(ij,ilev,var_local_n2n(nv3d+n))                                      

            pvec2d(ij,n,1:nbv)       = pvec3d(ij,ilev,var_local_n2n(nv3d+n),1:nbv)
            pmat2d(ij,n,1:nbv,1:nbv) = pmat3d(ij,ilev,var_local_n2n(nv3d+n),1:nbv,1:nbv) 
            peff2d(ij,n)             = peff3d(ij,ilev,var_local_n2n(nv3d+n))
          ELSE IF(var_local_n2n(nv3d+n) < nv3d+n) THEN
            trans(:,:,nv3d+n)        = trans(:,:,var_local_n2n(nv3d+n))
            work2d(ij,n)             = work2d(ij,var_local_n2n(nv3d+n)-nv3d)
            wvec2d(ij,n,1:nbv)       = wvec2d(ij,var_local_n2n(nv3d+n)-nv3d,1:nbv)        
            wmat2d(ij,n,1:nbv,1:nbv) = wmat2d(ij,var_local_n2n(nv3d+n)-nv3d,1:nbv,1:nbv)  
            nobs2d(ij,n)             = nobs2d(ij,var_local_n2n(nv3d+n)-nv3d)              
            
            pvec2d(ij,n,1:nbv)       = pvec2d(ij,var_local_n2n(nv3d+n)-nv3d,1:nbv)
            pmat2d(ij,n,1:nbv,1:nbv) = pmat2d(ij,var_local_n2n(nv3d+n)-nv3d,1:nbv,1:nbv)
            peff2d(ij,n)             = peff2d(ij,var_local_n2n(nv3d+n)-nv3d)
          ELSE
            CALL obs_local(ij,ilev,nv3d+n,hdxf,rdiag,rloc,dep,nobsl,logpfm)
            parm = work2d(ij,n)
            IF( dastype==3 ) parm = GAMMA_GMPF     !! LPFGM ; Walter and Potthast (2019)

            nobs2d(ij,n) = REAL( nobsl, r_size )
            CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,nv3d+n),logic_lpfgm, &
                            wvec2d(ij,n,1:nbv),wmat2d(ij,n,1:nbv,1:nbv)                             )
            IF( dastype==0 .or. dastype==1 ) work2d(ij,n) = parm !! LETKF, UPDATE

            IF( dastype==2 .or. dastype==3 ) parm = work2d(ij,n) ! LAPF or LPFGM
            CALL lpf_core  (nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,gr_common,type_pfmtx, &
                            asis2d(ij,n,1:nbv),peff2d(ij,n),pvec2d(ij,n,1:nbv),pmat2d(ij,n,1:nbv,1:nbv))
            IF( dastype==2 .or. dastype==3 ) work2d(ij,n) = parm ! LAPF or LPFGM
          END IF
        END DO
      END IF
    END DO !(End of IJ loop)   
    ptimer = MPI_WTIME()
    WRITE(6,'(A,2F10.2,A,i2)') '### TIMER in DAS_LETKF (MAIN LETKF):',ptimer-ptimer01,ptimer-ptimer00, "  Level: ",ilev
    ptimer01 = ptimer

!-----------------NOTE ---------------------------------------------------------------------------!
!##   wvec :: transform vector    ; LETKF for mean (IF logic_lpfgm=.TRUE. THEN wvec=0.0d0 )
!##   wmat :: transform matrix    ; LETKF for ptb  (IF logic_lpfgm=.TRUE. THEN wvec=transform matrix for dXb as LETKF)
!##   pvec :: weight of paticles  ; LPF   for mean
!##   pmat :: transform mtx       ; LPF   for particles  (!!CAUTION!! do not use pvec and pmat SIMULTANEOUSLY)
!##   asis :: SIS(log-likelihood) ; LPF   (inout)
!##   peff :: effective particle  ; LPF
!##   rmsp :: resampled or not    ; 0: no,  1: yes

!-----------------0. WEIGHT TREATMENTS       ---------------------------------------------------------------------------!
    ! currently no weight interpolation and weight smoother

!-----------------1. JUDGEMENT OF RESAMPLING ---------------------------------------------------------------------------!
    IF( dastype==1 .or. dastype==2 .or. dastype==3 )THEN
      !===> transform matrix = identity matrix I if no resampling
      DO ij=1,nij1      
          DO n=1,nv3d
            IF( peff3d(ij,ilev,n) >= resample_m  .or.  peff3d(ij,ilev,n) <= 0.0d0      )THEN !==> no  resampling              
              rsmp3d(ij,ilev,n)             = 0.0d0
              pmat3d(ij,ilev,n,1:nbv,1:nbv) = Ident(1:nbv,1:nbv)
            ELSE                                                                             !==> yes resampling
              asis3d(ij,ilev,n,1:nbv)       = 1.0d0/dble(nbv)
            ENDIF
          END DO
        IF( ilev==1 )THEN
          DO n=1,nv2d
            IF( peff2d(ij,     n) >= resample_m  .or.  peff2d(ij,     n) <= 0.0d0      )THEN
              !==> no  resampling
              rsmp2d(ij,     n)             = 0.0d0
              pmat2d(ij,     n,1:nbv,1:nbv) = Ident(1:nbv,1:nbv)
            ELSE
              !==> yes resampling
              asis2d(ij,     n,1:nbv)       = 1.0d0/dble(nbv)
            ENDIF
          END DO
        END IF
      END DO
    END IF

!-----------------2. LETKF, LPF, LPFGM CONTROLLER ------------------------------------------------------------------------------!
    DO ij=1,nij1  ! assumuming no variable loc.
      IF( dastype==0 )THEN      ! update :: LETKF for mean (wvec) & LETKF for ptb (wmat)
          asis3d(ij,ilev,1:nv3d,1:nbv) = 1.0d0/dble(nbv) ! no weight succession
          IF( ilev==1 ) &
          asis2d(ij,     1:nv2d,1:nbv) = 1.0d0/dble(nbv) ! no weight succession
      ELSE IF( dastype==1 )THEN ! update :: LPF   for mean (pvec) & LETKF for ptb (wmat)
          wvec3d(ij,ilev,1:nv3d,1:nbv) = pvec3d(ij,ilev,1:nv3d,1:nbv)
          asis3d(ij,ilev,1:nv3d,1:nbv) = 1.0d0/dble(nbv) ! no weight succession
          rsmp3d(ij,ilev,1:nv3d)       = 1.0d0           ! forced to resample
        IF( ilev==1 )THEN
          wvec2d(ij,     1:nv2d,1:nbv) = pvec2d(ij,     1:nv2d,1:nbv)
          asis2d(ij,     1:nv2d,1:nbv) = 1.0d0/dble(nbv) ! no weight succession
          rsmp2d(ij,     1:nv2d)       = 1.0d0           ! forced to resample
        ENDIF
      ELSE IF( dastype==2 )THEN ! update ::                       & LPF   for ensemble (pmat)
          wvec3d(ij,ilev,1:nv3d,1:nbv)       = 0.0d0  
          wmat3d(ij,ilev,1:nv3d,1:nbv,1:nbv) = pmat3d(ij,ilev,1:nv3d,1:nbv,1:nbv)
        IF( ilev==1 )THEN
          wvec2d(ij,     1:nv2d,1:nbv)       = 0.0d0
          wmat2d(ij,     1:nv2d,1:nbv,1:nbv) = pmat2d(ij,     1:nv2d,1:nbv,1:nbv)
        ENDIF
      ELSE IF ( dastype==3 )THEN ! update ::                        LPFGM for ensemble (wmat*pmat)
          !===> wvec is already defined to be 0.0 in common_letkf.f90
          DO n=1,nv3d
            wmat3d(ij,ilev,n,1:nbv,1:nbv) = matmul( wmat3d(ij,ilev,n,1:nbv,1:nbv),pmat3d(ij,ilev,n,1:nbv,1:nbv) )
          END DO
        IF( ilev==1 ) THEN
          DO n=1,nv2d
            wmat2d(ij,     n,1:nbv,1:nbv) = matmul( wmat2d(ij,     n,1:nbv,1:nbv),pmat2d(ij,     n,1:nbv,1:nbv) ) 
          ENDDO
        ENDIF
      END IF
    END DO
    ! MEMO :: LPF has only W matrix, no need to have w for mean updates

    !--------------3). ENSEMBLE UPDATES ------------------------------------------------------------------!
    DO ij=1,nij1
      DO n=1,nv3d
        DO m=1,nbv
          trans3(1:nbv,m)=wmat3d(ij,ilev,n,1:nbv,m)  + wvec3d(ij,ilev,n,1:nbv)
        END DO
        DO m=1,nbv
          anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n)
          DO k=1,nbv
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
                              & + gues3d(ij,ilev,k,n) * trans3(k,m)
          END DO
        END DO

        !==> RELAXATION (RTPP or RTPS)        
        IF( alph_relax/=0.0d0 .and. type_relax>=1) THEN
          tmpave        = sum( anal3d(ij,ilev,1:nbv,n) ) / dble( nbv )
          tmpptb(1:nbv) = anal3d(ij,ilev,1:nbv,n) - tmpave
          IF(      type_relax==1 ) THEN ! RTPP
            anal3d(ij,ilev,1:nbv,n) = tmpave + alph_relax *gues3d(ij,ilev,1:nbv,n) + &
                                     ( 1.0d0 - alph_relax)*tmpptb(1:nbv) 
          ELSE IF( type_relax==2 ) THEN ! RTPS
            gusspr = dsqrt( sum( gues3d(ij,ilev,1:nbv,n)**2.0d0 ) / REAL(nbv-1,r_size))
            anlspr = dsqrt( sum( tmpptb(        1:nbv  )**2.0d0 ) / REAL(nbv-1,r_size))
            IF( anlspr<gusspr .and. anlspr>0.0d0 )THEN
              tmpinf                  = 1.0d0 - alph_relax + alph_relax*(gusspr/anlspr) ! inflation factor
              anal3d(ij,ilev,1:nbv,n) = tmpave + tmpptb(1:nbv)*tmpinf
            END IF             
          END IF                                      
        ENDIF
      END DO
      IF(ilev >= 5) THEN !no analysis for upper-level Q
         DO m=1,nbv
            anal3d(ij,ilev,m,iv3d_q) = mean3d(ij,ilev,iv3d_q) &
                                   & + gues3d(ij,ilev,m,iv3d_q)
        END DO
      END IF

      IF(ilev == 1) THEN
        DO n=1,nv2d
          DO m=1,nbv
            trans2(1:nbv,m)=wmat2d(ij,n,1:nbv,m)  + wvec2d(ij,n,1:nbv)
          END DO
          DO m=1,nbv
            anal2d(ij,m,n)  = mean2d(ij,n)
            DO k=1,nbv
              anal2d(ij,m,n) = anal2d(ij,m,n) &
                           & + gues2d(ij,k,n) * trans2(k,m)
            END DO
          END DO

          !==> RELAXATION (RTPP or RTPS)
          IF( alph_relax/=0.0d0 .and. type_relax>=1 )THEN
            tmpave        = sum( anal2d(ij,1:nbv,n) ) / dble( nbv )
            tmpptb(1:nbv) = anal2d(ij,1:nbv,n) - tmpave
            IF(      type_relax==1 ) THEN ! RTPP
              anal2d(ij,1:nbv,n) = tmpave + alph_relax *gues2d(ij,1:nbv,n) + &
                                  ( 1.0d0 - alph_relax)*tmpptb(1:nbv) 
            ELSE IF( type_relax==2 ) THEN ! RTPS
              gusspr = dsqrt( sum( gues2d(ij,1:nbv,n)**2.0d0 ) / REAL(nbv-1,r_size))
              anlspr = dsqrt( sum( tmpptb(   1:nbv  )**2.0d0 ) / REAL(nbv-1,r_size))
              IF( anlspr<gusspr .and. anlspr>0.0d0 )THEN
                tmpinf             = 1.0d0 - alph_relax + alph_relax*(gusspr/anlspr) ! inflation factor
                anal2d(ij,1:nbv,n) = tmpave + tmpptb(1:nbv)*tmpinf
              END IF             
            END IF
          ENDIF
        END DO
      END IF
    END DO ! ij

    ptimer = MPI_WTIME()
    WRITE(6,'(A,2F10.2,A,i2)') '### TIMER in DAS_LETKF (WHT INTERP):',ptimer-ptimer01,ptimer-ptimer00, "  Level: ",ilev
    ptimer01 = ptimer
  !------------------ END OF NEW CODING FOR GATHER, FILTER, SCATTER, ANALYZE-----------------------------------------------!
  END DO   !(End of Assimilation loop; ilev)
  ptimer = MPI_WTIME()
  WRITE(6,'(A,2F10.2)') '### TIMER in DAS_LETKF (END LETKFs):',ptimer-ptimer01,ptimer-ptimer00
  ptimer01 = ptimer

  !!!!===> check consistency of random number
  !!!call com_rand_seed(nbv,iseed,srand)
  !!!print '(i5,5f8.5)', myrank, srand(1:5) 

  DEALLOCATE(hdxf,rdiag,rloc,dep)
  IF(cov_infl_mul < 0.0d0) THEN
    CALL gather_grd_mpi(0,work3d,work2d,work3dg,work2dg)
    IF(myrank == 0) THEN
      write(inflfile_o(12:21), '(I10.10)') ymdh
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. ',inflfile_o
      CALL write_grd4(inflfile_o,work3dg,work2dg)
    END IF
    DEALLOCATE(work3dg,work2dg,work3d,work2d)
  END IF
  !
  ! Effective Particle Size & Selectited Particle Numbers
  !
  ALLOCATE( work3dg(nlon,nlat,nlev,nv3d), work2dg(nlon,nlat,nv2d) )
  ALLOCATE( work3d(nij1,nlev,nv3d)      , work2d(nij1,nv2d)       )
    !==> effective particle size
    work3d = peff3d ; work2d = peff2d
    CALL gather_grd_mpi(0,work3d,work2d,work3dg,work2dg)
    IF(myrank == 0) THEN
      write(pefffile_o(12:21), '(I10.10)') ymdh
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. ',pefffile_o
      CALL write_grd4(pefffile_o,work3dg,work2dg)
    END IF

    !==> resampled or not (0: no, 1: yes)
    work3d = rsmp3d ; work2d = rsmp2d
    CALL gather_grd_mpi(0,work3d,work2d,work3dg,work2dg)
    IF(myrank == 0) THEN
      write(rsmpfile_o(12:21), '(I10.10)') ymdh
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. ',rsmpfile_o
      CALL write_grd4(rsmpfile_o,work3dg,work2dg)
    END IF
  DEALLOCATE(work3dg,work2dg,work3d,work2d)

  !
  ! Additive inflation
  !
  IF(sp_infl_add > 0.0d0) THEN
    CALL read_ens_mpi('addi',nbv,gues3d,gues2d)
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    CALL ensmean_grd(nbv,nij1,gues3d,gues2d,work3d,work2d)
    DO n=1,nv3d
      DO m=1,nbv
        DO k=1,nlev
          DO i=1,nij1
            gues3d(i,k,m,n) = gues3d(i,k,m,n) - work3d(i,k,n)
          END DO
        END DO
      END DO
    END DO
    DO n=1,nv2d
      DO m=1,nbv
        DO i=1,nij1
          gues2d(i,m,n) = gues2d(i,m,n) - work2d(i,n)
        END DO
      END DO
    END DO

    DEALLOCATE(work3d,work2d)
    WRITE(6,'(A)') '===== Additive covariance inflation ====='
    WRITE(6,'(A,F10.4)') '  parameter:',sp_infl_add
    WRITE(6,'(A)') '========================================='
!    parm = 0.7d0
!    DO ilev=1,nlev
!      parm_infl_damp(ilev) = 1.0d0 + parm &
!        & + parm * REAL(1-ilev,r_size)/REAL(nlev_dampinfl,r_size)
!      parm_infl_damp(ilev) = MAX(parm_infl_damp(ilev),1.0d0)
!    END DO
    DO n=1,nv3d
      DO m=1,nbv
        DO ilev=1,nlev
          DO ij=1,nij1
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
              & + gues3d(ij,ilev,m,n) * sp_infl_add
          END DO
        END DO
      END DO
    END DO
    DO n=1,nv2d
      DO m=1,nbv
        DO ij=1,nij1
          anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,m,n) * sp_infl_add
        END DO
      END DO
    END DO
  END IF

  !--> A.P 6/6/2018 - To save weight vector of LETKF
  IF( logic_wout )THEN
    ALLOCATE( work3d(nij1,nlev,nv3d), work3dg(nlon,nlat,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d)     , work2dg(nlon,nlat,nv2d)      )

    !==> max weight
    DO n=1,nv3d ; DO ilev=1,nlev ; DO ij=1,nij1
      work3d(ij,ilev,n) = maxval( abs(wvec3d(ij,ilev,n,:)) )
      !work3d(ij,ilev,n) =         sum(wmat3d(ij,ilev,n,:,1)) ! debug
    END DO      ; END DO         ; END DO
    DO n=1,nv2d ;                  DO ij=1,nij1
      work2d(ij,n)      = maxval( abs(wvec2d(ij,n,:)     ) )
    END DO      ;                  END DO
    CALL gather_grd_mpi(0,work3d,work2d,work3dg,work2dg)
    IF(myrank == 0) THEN
        WRITE(wvecfile,'(a7,i10.10,a)') './wvec/',ymdh,'_Mensmax.grd'
        WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. (max member)',wvecfile
        CALL write_grd4_1atm(wvecfile,work3dg,work2dg)
    END IF

    !==> each member for monitor (resulting weights)
     !debug!DO m=1,nbv
     DO m=1,2
      work3d(:,:,:) = wvec3d(:,:,:,m)
      work2d(:,:)   = wvec2d(:,:,m)

      !work3d(:,:,:) = wmat3d(:,:,:,m,1) !debug!
      CALL gather_grd_mpi(0,work3d,work2d,work3dg,work2dg)
       
      IF(myrank == 0) THEN
        WRITE(wvecfile,'(a7,i10.10,a2,i6.6,a4)') './wvec/',ymdh,'_M',m,'.grd'
        IF( m==1   ) WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. (fst member): ',wvecfile
        IF( m==nbv ) WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. (fnl member): ',wvecfile
        CALL write_grd4_1atm(wvecfile,work3dg,work2dg)
      END IF
    END DO
    DEALLOCATE(work3dg,work2dg,work3d,work2d)
  ENDIF

  !==> each member for SIS (need to save all member)
  IF( fgt_factor < 1.0d0 )THEN
    ALLOCATE( work3d(nij1,nlev,nv3d), work3dg(nlon,nlat,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d)     , work2dg(nlon,nlat,nv2d)      )
    DO m=1,nbv
      work3d(:,:,:) = asis3d(:,:,:,m)
      work2d(:,:)   = asis2d(:,:,m)
      CALL gather_grd_mpi(0,work3d,work2d,work3dg,work2dg)
       
      IF(myrank == 0) THEN
        WRITE(wvecfile,'(a7,i10.10,a2,i6.6,a4)') './asis_lpf/',ymdh,'_M',m,'.grd'
        IF( m==1   ) WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. (fst member): ',wvecfile
        IF( m==nbv ) WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. (fnl member): ',wvecfile
        CALL write_grd4_1atm(wvecfile,work3dg,work2dg)
      END IF
    END DO
    DEALLOCATE(work3dg,work2dg,work3d,work2d)
  ENDIF
  !--> END OF EDITS
!
  ptimer = MPI_WTIME()
  WRITE(6,'(A,2F10.2)') '### TIMER in DAS_LETKF (EXIT LETKF):',ptimer-ptimer01,ptimer-ptimer00
  ptimer01 = ptimer
!
  DEALLOCATE(nobs3d, nobs2d )
  DEALLOCATE(wvec3d, wvec2d )
  DEALLOCATE(wmat3d, wmat2d )
  DEALLOCATE(pvec3d, pvec2d )
  DEALLOCATE(pmat3d, pmat2d )
  DEALLOCATE(logpfm,mean3d,mean2d)
  RETURN
END SUBROUTINE das_letkf

!-----------------------------------------------------------------------
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
!-----------------------------------------------------------------------
SUBROUTINE obs_local(ij,ilev,nvar,hdxf,rdiag,rloc,dep,nobsl,logpfm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ij,ilev,nvar
  REAL(r_size),INTENT(IN) :: logpfm(nij1,nlev)
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT) :: nobsl
  REAL(r_size) :: minlon,maxlon,minlat,maxlat,dist,dlev
  REAL(r_size) :: tmplon,tmplat,tmperr,tmpwgt(nlev)
  INTEGER :: tmpqc
  INTEGER,ALLOCATABLE:: nobs_use(:)
!TVS  INTEGER,ALLOCATABLE:: ntvs_use_prof(:),ntvs_use_inst(:),ntvs_use_slot(:)
  INTEGER :: imin,imax,jmin,jmax,im,ichan
  INTEGER :: n,nn,tvnn,iobs
!
! INITIALIZE
!
  IF( nobs > 0 ) THEN
    ALLOCATE(nobs_use(nobs))
  END IF
!TVS  IF( ntvs > 0 ) THEN
!TVS    ALLOCATE(ntvs_use_prof(ntvs))
!TVS    ALLOCATE(ntvs_use_inst(ntvs))
!TVS    ALLOCATE(ntvs_use_slot(ntvs))
!TVS  END IF
!
! data search
!
  minlon = lon1(ij) - dlon_zero(ij)
  maxlon = lon1(ij) + dlon_zero(ij)
  minlat = lat1(ij) - dlat_zero
  maxlat = lat1(ij) + dlat_zero
  IF(maxlon - minlon >= 360.0d0) THEN
    minlon = 0.0d0
    maxlon = 360.0d0
  END IF

  DO jmin=1,nlat-2
    IF(minlat < lat(jmin+1)) EXIT
  END DO
  DO jmax=1,nlat-2
    IF(maxlat < lat(jmax+1)) EXIT
  END DO
  nn = 1
!TVS  tvnn = 1
  IF(minlon >= 0 .AND. maxlon <= 360.0) THEN
    DO imin=1,nlon-1
      IF(minlon < lon(imin+1)) EXIT
    END DO
    DO imax=1,nlon-1
      IF(maxlon < lon(imax+1)) EXIT
    END DO
    IF( nobs > 0 ) &
    & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS    IF( ntvs > 0 ) &
!TVS    & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS    &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
  ELSE IF(minlon >= 0 .AND. maxlon > 360.0) THEN
    DO imin=1,nlon-1
      IF(minlon < lon(imin+1)) EXIT
    END DO
    maxlon = maxlon - 360.0d0
    IF(maxlon > 360.0d0) THEN
      imin = 1
      imax = nlon
      IF( nobs > 0 ) &
      & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS      IF( ntvs > 0 ) &
!TVS      & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS      &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
    ELSE
      DO imax=1,nlon-1
        IF(maxlon < lon(imax+1)) EXIT
      END DO
      IF(imax > imin) THEN
        imin = 1
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      ELSE
        imin = 1
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
        DO imin=1,nlon-1
          IF(minlon < lon(imin+1)) EXIT
        END DO
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      END IF
    END IF
  ELSE IF(minlon < 0 .AND. maxlon <= 360.0d0) THEN
    DO imax=1,nlon-1
      IF(maxlon < lon(imax+1)) EXIT
    END DO
    minlon = minlon + 360.0d0
    IF(minlon < 0) THEN
      imin = 1
      imax = nlon
      IF( nobs > 0 ) &
      & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS      IF( ntvs > 0 ) &
!TVS      & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS      &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
    ELSE
      DO imin=1,nlon-1
        IF(minlon < lon(imin+1)) EXIT
      END DO
      IF(imin < imax) THEN
        imin = 1
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      ELSE
        imin = 1
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
        DO imin=1,nlon-1
          IF(minlon < lon(imin+1)) EXIT
        END DO
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      END IF
    END IF
  ELSE
    maxlon = maxlon - 360.0d0
    minlon = minlon + 360.0d0
    IF(maxlon > 360.0 .OR. minlon < 0) THEN
      imin = 1
      imax = nlon
      IF( nobs > 0 ) &
      & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS      IF( ntvs > 0 ) &
!TVS      & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS      &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
    ELSE
      DO imin=1,nlon-1
        IF(minlon < lon(imin+1)) EXIT
      END DO
      DO imax=1,nlon-1
        IF(maxlon < lon(imax+1)) EXIT
      END DO
      IF(imin > imax) THEN
        imin = 1
        imax = nlon
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      ELSE
        IF( nobs > 0 ) &
        & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!TVS        IF( ntvs > 0 ) &
!TVS        & CALL tvs_local_sub(imin,imax,jmin,jmax,tvnn, &
!TVS        &                    ntvs_use_prof,ntvs_use_inst,ntvs_use_slot)
      END IF
    END IF
  END IF
  nn = nn-1
!TVS  tvnn = tvnn -1
!TVS  IF( nn < 1 .AND. tvnn < 1 ) THEN
  IF(nn < 1) THEN
    nobsl = 0
    RETURN
  END IF
!
! CONVENTIONAL
!
  nobsl = 0
  IF(nn > 0) THEN
    DO n=1,nn
      !
      ! vertical localization
      !
      IF(NINT(obselm(nobs_use(n))) == id_ps_obs .AND. ilev > 1) THEN
        dlev = ABS(LOG(obsdat(nobs_use(n))) - logpfm(ij,ilev))
      ELSE IF(NINT(obselm(nobs_use(n))) /= id_ps_obs) THEN
        dlev = ABS(LOG(obslev(nobs_use(n))) - logpfm(ij,ilev))
      ELSE
        dlev = 0.0d0
      END IF
      IF(dlev > dist_zerov) CYCLE
      !
      ! horizontal localization
      !
      tmplon=obslon(nobs_use(n))
      tmplat=obslat(nobs_use(n))
      CALL com_distll_1( tmplon, tmplat,lon1(ij), lat1(ij), dist)
      IF(dist > dist_zero ) CYCLE
      !
      ! variable localization
      !
      SELECT CASE(NINT(obselm(nobs_use(n))))
      CASE(id_u_obs)
        iobs=1
      CASE(id_v_obs)
        iobs=2
      CASE(id_t_obs)
        iobs=3
      CASE(id_q_obs)
        iobs=4
      CASE(id_rh_obs)
        iobs=5
      CASE(id_ps_obs)
        iobs=6
      CASE(id_rain_obs)
        iobs=7
      END SELECT
      IF(var_local(nvar,iobs) < TINY(var_local)) CYCLE

      nobsl = nobsl + 1
      hdxf(nobsl,:) = obshdxf(nobs_use(n),:)
      dep(nobsl)    = obsdep(nobs_use(n))
      !
      ! Observational localization
      !
      tmperr=obserr(nobs_use(n))
      rdiag(nobsl) = tmperr * tmperr
      rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) &
                  & * var_local(nvar,iobs)
    END DO
  END IF
!TVS!
!TVS! ATOVS
!TVS!
!TVS  IF(tvnn > 0) THEN
!TVS    DO n=1,tvnn
!TVS      tmplon=tvslon(ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))
!TVS      tmplat=tvslat(ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))
!TVS      CALL com_distll_1( tmplon, tmplat, lon1(ij), lat1(ij), dist)
!TVS      IF( dist > dist_zero) CYCLE
!TVS
!TVS      DO ichan=1,ntvsch(ntvs_use_inst(n))
!TVS        tmperr=tvserr(ichan,ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))
!TVS        tmpqc=tvsqc(ichan,ntvs_use_prof(n),ntvs_use_inst(n),ntvs_use_slot(n))
!TVS        tmpwgt(:)=tvswgt(:,ichan, &
!TVS                         & ntvs_use_prof(n), &
!TVS                         & ntvs_use_inst(n), &
!TVS                         & ntvs_use_slot(n))
!TVS        IF( tmpqc == 1 .AND. tmpwgt(ilev) > 0.05D0 ) THEN
!TVS          nobsl = nobsl + 1
!TVS          DO im = 1, nbv
!TVS            hdxf(nobsl,im) = tvshdxf(im,ichan, &
!TVS                              & ntvs_use_prof(n), &
!TVS                              & ntvs_use_inst(n), &
!TVS                              & ntvs_use_slot(n))
!TVS          END DO
!TVS          dep(nobsl)    = tvsdep(ichan, &
!TVS                              & ntvs_use_prof(n), &
!TVS                              & ntvs_use_inst(n), &
!TVS                              & ntvs_use_slot(n))
!TVS          rdiag(nobsl)  = tmperr * tmperr &
!TVS                        & * exp(0.5d0 * (dist/sigma_obs)**2) &
!TVS                        & / (tmpwgt(ilev) * tmpwgt(ilev))
!TVS        END IF
!TVS      END DO
!TVS    END DO
!TVS  END IF
!
! DEBUG
! IF( ILEV == 1 .AND. ILON == 1 ) &
! & WRITE(6,*) 'ILEV,ILON,ILAT,NN,TVNN,NOBSL=',ilev,ij,nn,tvnn,nobsl
!
  IF( nobsl > nobstotal ) THEN
    WRITE(6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL=',nobsl,' > NOBSTOTAL=',nobstotal
    WRITE(6,*) 'IJ,NN,TVNN=', ij, nn, tvnn
    STOP 99
  END IF
!
  IF( nobs > 0 ) THEN
    DEALLOCATE(nobs_use)
  END IF
!TVS  IF( ntvs > 0 ) THEN
!TVS    DEALLOCATE(ntvs_use_prof)
!TVS    DEALLOCATE(ntvs_use_inst)
!TVS    DEALLOCATE(ntvs_use_slot)
!TVS  END IF
!
  RETURN
END SUBROUTINE obs_local

SUBROUTINE obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
  INTEGER,INTENT(INOUT) :: nn, nobs_use(nobs)
  INTEGER :: j,n,ib,ie,ip

  DO j=jmin,jmax
    IF(imin > 1) THEN
      ib = nobsgrd(imin-1,j)+1
    ELSE
      IF(j > 1) THEN
        ib = nobsgrd(nlon,j-1)+1
      ELSE
        ib = 1
      END IF
    END IF
    ie = nobsgrd(imax,j)
    n = ie - ib + 1
    IF(n == 0) CYCLE
    DO ip=ib,ie
      IF(nn > nobs) THEN
        WRITE(6,*) 'FATALERROR, NN > NOBS', NN, NOBS
      END IF
      nobs_use(nn) = ip
      nn = nn + 1
    END DO
  END DO

  RETURN
END SUBROUTINE obs_local_sub

!TVSSUBROUTINE tvs_local_sub(imin,imax,jmin,jmax,nn,ntvs_prof,ntvs_inst,ntvs_slot)
!TVS  INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
!TVS  INTEGER,INTENT(INOUT) :: nn, ntvs_prof(ntvs), ntvs_inst(ntvs), ntvs_slot(ntvs)
!TVS  INTEGER :: j,n,ib,ie,ip
!TVS  INTEGER :: islot, iinst
!TVS
!TVS  DO j=jmin,jmax
!TVS    DO islot=1,nslots
!TVS      DO iinst=1,ninstrument
!TVS        IF(imin > 1) THEN
!TVS          ib = ntvsgrd(imin-1,j,iinst,islot)+1
!TVS        ELSE
!TVS          IF(j > 1) THEN
!TVS            ib = ntvsgrd(nlon,j-1,iinst,islot)+1
!TVS          ELSE
!TVS            ib = 1
!TVS          END IF
!TVS        END IF
!TVS        ie = ntvsgrd(imax,j,iinst,islot)
!TVS        n = ie - ib + 1
!TVS        IF(n == 0) CYCLE
!TVS        DO ip=ib,ie
!TVS          IF(nn > nobs) THEN
!TVS            WRITE(6,*) 'FATALERROR, NN > NTVS', NN, NTVS
!TVS          END IF
!TVS          ntvs_prof(nn)=ip
!TVS          ntvs_inst(nn)=iinst
!TVS          ntvs_slot(nn)=islot
!TVS          nn = nn + 1
!TVS        END DO
!TVS      END DO
!TVS    END DO
!TVS  END DO
!TVS  RETURN
!TVSEND SUBROUTINE tvs_local_sub
!TVS!-----------------------------------------------------------------------
!TVS! Data Assimilation for VARBC
!TVS!-----------------------------------------------------------------------
!TVSSUBROUTINE das_vbc(um,vm,tm,qm,qlm,psm,vbcf,vbca)
!TVS  USE common_mtx
!TVS  IMPLICIT NONE
!TVS  REAL(r_size),INTENT(IN) :: um(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: vm(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: tm(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: qm(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: qlm(nij1,nlev)
!TVS  REAL(r_size),INTENT(IN) :: psm(nij1)
!TVS  REAL(r_size),INTENT(INOUT) :: vbcf(maxvbc,maxtvsch,ninstrument)
!TVS  REAL(r_size),INTENT(OUT)   :: vbca(maxvbc,maxtvsch,ninstrument)
!TVS  REAL(r_sngl) :: u4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: v4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: t4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: q4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: ql4(nlon,nlat,nlev)
!TVS  REAL(r_sngl) :: ps4(nlon,nlat)
!TVS  REAL(r_size) :: u(nlon,nlat,nlev)
!TVS  REAL(r_size) :: v(nlon,nlat,nlev)
!TVS  REAL(r_size) :: t(nlon,nlat,nlev)
!TVS  REAL(r_size) :: q(nlon,nlat,nlev)
!TVS  REAL(r_size) :: ql(nlon,nlat,nlev)
!TVS  REAL(r_size) :: ps(nlon,nlat)
!TVS  REAL(r_size) :: p_full(nlon,nlat,nlev)
!TVS  REAL(r_size),ALLOCATABLE :: hx(:,:,:,:)
!TVS  REAL(r_size),ALLOCATABLE :: pred(:,:,:,:,:)
!TVS  INTEGER,ALLOCATABLE :: tmpqc(:,:,:)
!TVS  REAL(r_size),ALLOCATABLE :: tmpwgt(:,:,:,:)
!TVS  REAL(r_size) :: a(maxvbc,maxvbc)
!TVS  REAL(r_size) :: b(maxvbc)
!TVS  REAL(r_size) :: ainv(maxvbc,maxvbc)
!TVS  INTEGER:: ntvschan1(maxtvsch,ninstrument)
!TVS  INTEGER:: i,j,k,n,islot,nn
!TVS
!TVS  PRINT *,'Hello from das_vbc'
!TVS
!TVS  IF(ntvs == 0) THEN
!TVS    PRINT *,'No radiance data: das_vbc skipped..'
!TVS!$OMP PARALLEL WORKSHARE
!TVS    vbca = vbcf
!TVS!$OMP END PARALLEL WORKSHARE
!TVS    RETURN
!TVS  END IF
!TVS
!TVS  CALL gather_grd_mpi(0,um,vm,tm,qm,qlm,psm,u4,v4,t4,q4,ql4,ps4)
!TVS  n = nlon*nlat*nlev
!TVS  CALL MPI_BARRIER(MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(u4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(v4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(t4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(q4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  CALL MPI_BCAST(ql4(1,1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS  n = nlon*nlat
!TVS  CALL MPI_BCAST(ps4(1,1),n,MPI_REAL,0,MPI_COMM_WORLD,i)
!TVS!$OMP PARALLEL WORKSHARE
!TVS  u = REAL(u4,r_size)
!TVS  v = REAL(v4,r_size)
!TVS  t = REAL(t4,r_size)
!TVS  q = REAL(q4,r_size)
!TVS  ql = REAL(ql4,r_size)
!TVS  ps = REAL(ps4,r_size)
!TVS!$OMP END PARALLEL WORKSHARE
!TVS  CALL calc_pfull(ps,p_full)
!TVS
!TVS  ALLOCATE( hx(maxtvsch,maxtvsprof,ninstrument,nslots) )
!TVS  ALLOCATE( pred(maxvbc,maxtvsch,maxtvsprof,ninstrument,nslots) )
!TVS  ALLOCATE( tmpqc(maxtvsch,maxtvsprof,ninstrument) )
!TVS  ALLOCATE( tmpwgt(nlev,maxtvsch,maxtvsprof,ninstrument) )
!TVS  DO islot=1,nslots
!TVS!    IF(SUM(ntvsprofslots(:,islot)) == 0) CYCLE
!TVS    ntvsprof(:) = ntvsprofslots(:,islot)
!TVS    CALL Trans_XtoY_tvs(u,v,t,q,ql,ps,p_full, &
!TVS      & tvslon(:,:,islot),tvslat(:,:,islot),tvszenith(:,:,islot),&
!TVS      & tvsskin(:,:,islot),tvsstmp(:,:,islot),tvsclw(:,:,islot),&
!TVS      & tvsemis(:,:,:,islot),tmpqc,hx(:,:,:,islot),tmpwgt,pred(:,:,:,:,islot))
!TVS  END DO
!TVS  DEALLOCATE(tmpqc,tmpwgt)
!TVS
!TVS!$OMP PARALLEL PRIVATE(j,k,n,a,b,ainv)
!TVS!$OMP WORKSHARE
!TVS  vbca = 0.0d0
!TVS!$OMP END WORKSHARE
!TVS!$OMP DO SCHEDULE(DYNAMIC)
!TVS  DO k=1,ninstrument
!TVS    DO j=1,maxtvsch
!TVS      !
!TVS      ! Parallel processing
!TVS      !
!TVS      IF(MOD(j+maxtvsch*(k-1)-1,nprocs) /= myrank) CYCLE
!TVS      !
!TVS      ! DATA NUMBER
!TVS      !
!TVS      ntvschan(j,k) = SUM(tvsqc(j,:,k,:))
!TVS      IF(msw_vbc .AND. ntvschan(j,k) /= 0 ) THEN
!TVS        PRINT '(3A,I3,A,I6)',' >> VBC executed for instrument,channel,ntvsl: ',&
!TVS                            & tvsname(k),',',tvsch(j,k),',',ntvschan(j,k)
!TVS        CALL vbc_local(j,k,ntvschan(j,k),hx,pred,a,b)
!TVS        CALL mtx_inv(maxvbc,a,ainv)
!TVS        vbca(:,j,k) = vbcf(:,j,k)
!TVS        DO n=1,maxvbc
!TVS          vbca(:,j,k) = vbca(:,j,k) - ainv(:,n)*b(n) !ATTN: sign for beta
!TVS        END DO
!TVS      ELSE
!TVS        PRINT '(3A,I3,A,I6)',' !! NO VBC executed for instrument,channel,ntvsl: ',&
!TVS                            & tvsname(k),',',tvsch(j,k),',',ntvschan(j,k)
!TVS        vbca(:,j,k) = vbcf(:,j,k)
!TVS      END IF
!TVS    END DO
!TVS  END DO
!TVS!$OMP END DO
!TVS!$OMP WORKSHARE
!TVS  vbcf = vbca
!TVS  ntvschan1 = ntvschan
!TVS!$OMP END WORKSHARE
!TVS!$OMP END PARALLEL
!TVS  DEALLOCATE(hx,pred)
!TVS  n = maxvbc*maxtvsch*ninstrument
!TVS  CALL MPI_BARRIER(MPI_COMM_WORLD,j)
!TVS  CALL MPI_ALLREDUCE(vbcf,vbca,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,j)
!TVS  n = maxtvsch*ninstrument
!TVS  CALL MPI_BARRIER(MPI_COMM_WORLD,j)
!TVS  CALL MPI_ALLREDUCE(ntvschan1,ntvschan,n,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,j)
!TVS
!TVS  RETURN
!TVSEND SUBROUTINE das_vbc
!TVS!-----------------------------------------------------------------------
!TVS!  (in) ichan: channnel
!TVS!  (in) iinst: sensor
!TVS!  (out) a = B_beta^-1 + p R^-1 p^T
!TVS!  (out) b = p R^-1 d
!TVS!-----------------------------------------------------------------------
!TVSSUBROUTINE vbc_local(ichan,iinst,ntvsl,hx,pred,a,b)
!TVS  IMPLICIT NONE
!TVS  INTEGER,PARAMETER :: msw=1
!TVS  INTEGER,PARAMETER :: nmin=400
!TVS  INTEGER,INTENT(IN) :: ichan,iinst,ntvsl
!TVS  REAL(r_size),INTENT(IN) :: hx(maxtvsch,maxtvsprof,ninstrument,nslots)
!TVS  REAL(r_size),INTENT(IN) :: pred(maxvbc,maxtvsch,maxtvsprof,ninstrument,nslots)
!TVS  REAL(r_size),INTENT(OUT) :: a(maxvbc,maxvbc)
!TVS  REAL(r_size),INTENT(OUT) :: b(maxvbc)
!TVS  REAL(r_size) :: dep,dep0
!TVS  REAL(r_size) :: bias,bias0
!TVS  REAL(r_size) :: r,tmp
!TVS  INTEGER:: islot, iprof, i,j,n
!TVS
!TVS  a = 0.0d0
!TVS  b = 0.0d0
!TVS  dep = 0.0d0
!TVS  dep0 = 0.0d0
!TVS  bias = 0.0d0
!TVS  bias0 = 0.0d0
!TVS  n = 0
!TVS  DO islot=1,nslots
!TVS    DO iprof=1,maxtvsprof
!TVS      IF(tvsqc(ichan,iprof,iinst,islot)/=1) CYCLE
!TVS      !
!TVS      ! R
!TVS      !
!TVS      r = tvserr(ichan,iprof,iinst,islot)**2
!TVS      !
!TVS      ! p R^-1 p^T
!TVS      !
!TVS      DO j=1,maxvbc
!TVS        DO i=1,maxvbc
!TVS          a(i,j) = a(i,j) &
!TVS               & + pred(i,ichan,iprof,iinst,islot) &
!TVS               & * pred(j,ichan,iprof,iinst,islot) / r
!TVS        END DO
!TVS      END DO
!TVS      !
!TVS      ! B_beta^-1
!TVS      !
!TVS      IF(msw == 1) THEN ! Y.Sato
!TVS        IF(ntvsl < nmin) THEN
!TVS          tmp = REAL(nmin,r_size) / r
!TVS
!TVS        ELSE
!TVS          tmp = (REAL(ntvsl,r_size) &
!TVS            & / (LOG10(REAL(ntvsl,r_size)/REAL(nmin,r_size))+1.0d0)) / r
!TVS        END IF
!TVS      ELSE IF(msw == 2) THEN ! D.Dee
!TVS        tmp = REAL(ntvsl,r_size) / r
!TVS      ELSE ! Constant
!TVS        tmp = 100.0d0
!TVS      END IF
!TVS      DO i=1,maxvbc
!TVS        a(i,i) = a(i,i) + tmp
!TVS      END DO
!TVS      !
!TVS      ! p R^-1 d
!TVS      !
!TVS      b(:) = b(:) + pred(:,ichan,iprof,iinst,islot) / r &
!TVS                & *(tvsdat(ichan,iprof,iinst,islot)-hx(ichan,iprof,iinst,islot))
!TVS      bias = bias+tvsdat(ichan,iprof,iinst,islot)-hx(ichan,iprof,iinst,islot)
!TVS      dep = dep+(tvsdat(ichan,iprof,iinst,islot)-hx(ichan,iprof,iinst,islot))**2
!TVS      bias0= bias0+tvsdep(ichan,iprof,iinst,islot)
!TVS      dep0= dep0+tvsdep(ichan,iprof,iinst,islot)**2
!TVS      n = n+1
!TVS    END DO
!TVS  END DO
!TVS
!TVS  dep = SQRT(dep / REAL(n,r_size))
!TVS  dep0 = SQRT(dep0 / REAL(n,r_size))
!TVS  bias = bias / REAL(n,r_size)
!TVS  bias0 = bias0 / REAL(n,r_size)
!TVS  PRINT '(2A,I3,4F12.4)',' >> D monit: ',tvsname(iinst),tvsch(ichan,iinst),bias0,bias,dep0,dep
!TVS
!TVS  RETURN
!TVSEND SUBROUTINE vbc_local

END MODULE letkf_tools
