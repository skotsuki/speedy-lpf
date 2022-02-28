MODULE letkf_obs
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   04/03/2013 Takemasa MIYOSHI  separating obs operator
!   01/23/2009 Takemasa MIYOSHI  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_speedy
  USE common_obs_speedy
  USE common_mpi_speedy
  USE common_letkf

  IMPLICIT NONE
  PUBLIC

  INTEGER,SAVE :: nobs
  INTEGER,PARAMETER :: nslots=1 ! number of time slots for 4D-LETKF
  INTEGER,PARAMETER :: nbslot=1 ! basetime slot
!  REAL(r_size),PARAMETER :: sigma_obs=500.0d3
!  REAL(r_size),PARAMETER :: sigma_obsv=0.4d0
  REAL(r_size),PARAMETER :: sigma_obsv=0.1d0
  REAL(r_size),PARAMETER :: sigma_obst=3.0d0
  REAL(r_size),SAVE :: dist_zero
  REAL(r_size),SAVE :: dist_zerov
  REAL(r_size),ALLOCATABLE,SAVE :: dlon_zero(:)
  REAL(r_size),SAVE :: dlat_zero
  REAL(r_size),ALLOCATABLE,SAVE :: obselm(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslon(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslev(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obserr(:)
!  REAL(r_size),ALLOCATABLE,SAVE :: obsk(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdep(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obshdxf(:,:)
  INTEGER,SAVE :: nobsgrd(nlon,nlat)

!  REAL(r_size),SAVE :: sigma_obs=9999.0d3  

CONTAINS
!-----------------------------------------------------------------------
! Initialize
!-----------------------------------------------------------------------
SUBROUTINE set_letkf_obs
  IMPLICIT NONE
  REAL(r_size),PARAMETER :: gross_error=10.0d0
  REAL(r_size) :: dlon1,dlon2,dlon,dlat
  REAL(r_size),ALLOCATABLE :: wk2d(:,:)
  INTEGER,ALLOCATABLE :: iwk2d(:,:)
  REAL(r_size),ALLOCATABLE :: tmpelm(:)
  REAL(r_size),ALLOCATABLE :: tmplon(:)
  REAL(r_size),ALLOCATABLE :: tmplat(:)
  REAL(r_size),ALLOCATABLE :: tmplev(:)
  REAL(r_size),ALLOCATABLE :: tmpdat(:)
  REAL(r_size),ALLOCATABLE :: tmperr(:)
  REAL(r_size),ALLOCATABLE :: tmpk(:)
  REAL(r_size),ALLOCATABLE :: tmpdep(:)
  REAL(r_size),ALLOCATABLE :: tmphdxf(:,:)
  INTEGER,ALLOCATABLE :: tmpqc0(:,:)
  INTEGER,ALLOCATABLE :: tmpqc(:)
  REAL(r_size),ALLOCATABLE :: tmp2elm(:)
  REAL(r_size),ALLOCATABLE :: tmp2lon(:)
  REAL(r_size),ALLOCATABLE :: tmp2lat(:)
  REAL(r_size),ALLOCATABLE :: tmp2lev(:)
  REAL(r_size),ALLOCATABLE :: tmp2dat(:)
  REAL(r_size),ALLOCATABLE :: tmp2err(:)
!  REAL(r_size),ALLOCATABLE :: tmp2k(:)
  REAL(r_size),ALLOCATABLE :: tmp2dep(:)
  REAL(r_size),ALLOCATABLE :: tmp2hdxf(:,:)
  INTEGER :: nobslots(nslots)
  INTEGER :: n,i,j,ierr,islot,nn,l,im
  INTEGER :: nj(0:nlat-1)
  INTEGER :: njs(1:nlat-1)
  CHARACTER(15) :: obsfile='obsTTNNNNNN.dat'

  WRITE(6,'(A)') 'Hello from set_letkf_obs'

  dist_zero = sigma_obs * SQRT(10.0d0/3.0d0) * 2.0d0
  dist_zerov = sigma_obsv * SQRT(10.0d0/3.0d0) * 2.0d0
  dlat_zero = dist_zero / pi / re * 180.0d0
  ALLOCATE(dlon_zero(nij1))
  DO i=1,nij1
    !debug by A. Yamazaki SOLA 2017!dlon_zero(i) = dlat_zero / COS(pi*lat1(i)/180.0d0)
    CALL search_longestlat(dist_zero,lat1(i),dlat_zero,dlon_zero(i))
    ! [in]  :: dist_zero, lat1(i), dlat_zero
    ! [out] :: dlon_zero(i)
  END DO

  IF(myrank == 0) THEN !Assuming all members have the identical obs records
    DO islot=1,nslots
      im = myrank+1
      WRITE(obsfile(4:11),'(I2.2,I6.6)') islot,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',obsfile
      CALL get_nobs(obsfile,8,nobslots(islot))
    END DO
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(nobslots,nslots,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nobs = SUM(nobslots)
  WRITE(6,'(I10,A)') nobs,' TOTAL OBSERVATIONS INPUT'

  IF(nobs == 0) THEN
    WRITE(6,'(A)') 'No observation assimilated'
    RETURN
  END IF
!
! INITIALIZE GLOBAL VARIABLES
!
  ALLOCATE( tmpelm(nobs) )
  ALLOCATE( tmplon(nobs) )
  ALLOCATE( tmplat(nobs) )
  ALLOCATE( tmplev(nobs) )
  ALLOCATE( tmpdat(nobs) )
  ALLOCATE( tmperr(nobs) )
  ALLOCATE( tmpk(nobs) )
  ALLOCATE( tmpdep(nobs) )
  ALLOCATE( tmphdxf(nobs,nbv) )
  ALLOCATE( tmpqc0(nobs,nbv) )
  ALLOCATE( tmpqc(nobs) )
  tmpqc0 = 0
  tmphdxf = 0.0d0
!
! reading observation data
!
  nn=0
  timeslots: DO islot=1,nslots
    IF(nobslots(islot) == 0) CYCLE
    l=0
    DO
      im = myrank+1 + nprocs * l
      IF(im > nbv) EXIT
      WRITE(obsfile(4:11),'(I2.2,I6.6)') islot,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',obsfile
      CALL read_obs2(obsfile,nobslots(islot),&
       & tmpelm(nn+1:nn+nobslots(islot)),tmplon(nn+1:nn+nobslots(islot)),&
       & tmplat(nn+1:nn+nobslots(islot)),tmplev(nn+1:nn+nobslots(islot)),&
       & tmpdat(nn+1:nn+nobslots(islot)),tmperr(nn+1:nn+nobslots(islot)),&
       & tmphdxf(nn+1:nn+nobslots(islot),im),tmpqc0(nn+1:nn+nobslots(islot),im))
      l = l+1
    END DO
    nn = nn + nobslots(islot)
  END DO timeslots

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ALLOCATE(wk2d(nobs,nbv))
  wk2d = tmphdxf
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(wk2d,tmphdxf,nobs*nbv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  DEALLOCATE(wk2d)
  ALLOCATE(iwk2d(nobs,nbv))
  iwk2d = tmpqc0
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(iwk2d,tmpqc0,nobs*nbv,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  DEALLOCATE(iwk2d)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i)
  DO n=1,nobs
    tmpqc(n) = MINVAL(tmpqc0(n,:))
    IF(tmpqc(n) /= 1) CYCLE
    tmpdep(n) = tmphdxf(n,1)
    DO i=2,nbv
      tmpdep(n) = tmpdep(n) + tmphdxf(n,i)
    END DO
    tmpdep(n) = tmpdep(n) / REAL(nbv,r_size)
    DO i=1,nbv
      tmphdxf(n,i) = tmphdxf(n,i) - tmpdep(n) ! Hdx
    END DO
    tmpdep(n) = tmpdat(n) - tmpdep(n) ! y-Hx
    IF(ABS(tmpdep(n)) > gross_error*tmperr(n)) THEN !gross error
      tmpqc(n) = 0
    END IF
  END DO
!$OMP END PARALLEL DO
  DEALLOCATE(tmpqc0)

  WRITE(6,'(I10,A)') SUM(tmpqc),' OBSERVATIONS TO BE ASSIMILATED'

  CALL monit_dep(nobs,tmpelm,tmpdep,tmpqc)
!
! temporal observation localization
!
!  nn = 0
!  DO islot=1,nslots
!    tmperr(nn+1:nn+nobslots(islot)) = tmperr(nn+1:nn+nobslots(islot)) &
!      & * exp(0.25d0 * (REAL(islot-nbslot,r_size) / sigma_obst)**2)
!    nn = nn + nobslots(islot)
!  END DO
!
! SELECT OBS IN THE NODE
!
  nn = 0
  DO n=1,nobs
    IF(tmpqc(n) /= 1) CYCLE
!    IF(tmplat(n) < MINVAL(lat1) .OR. MAXVAL(lat1) < tmplat(n)) THEN
!      dlat = MIN( ABS(MINVAL(lat1)-tmplat(n)),ABS(MAXVAL(lat1)-tmplat(n)) )
!      IF(dlat > dlat_zero) CYCLE
!    END IF
!    IF(tmplon(n) < MINVAL(lon1) .OR. MAXVAL(lon1) < tmplon(n)) THEN
!      dlon1 = ABS(MINVAL(lon1) - tmplon(n))
!      dlon1 = MIN(dlon1,360.0d0-dlon1)
!      dlon2 = ABS(MAXVAL(lon1) - tmplon(n))
!      dlon2 = MIN(dlon2,360.0d0-dlon2)
!      dlon =  MIN(dlon1,dlon2) &
!         & * pi*re*COS(tmplat(n)*pi/180.d0)/180.0d0
!      IF(dlon > dist_zero) CYCLE
!    END IF
    nn = nn+1
    tmpelm(nn) = tmpelm(n)
    tmplon(nn) = tmplon(n)
    tmplat(nn) = tmplat(n)
    tmplev(nn) = tmplev(n)
    tmpdat(nn) = tmpdat(n)
    tmperr(nn) = tmperr(n)
    tmpk(nn) = tmpk(n)
    tmpdep(nn) = tmpdep(n)
    tmphdxf(nn,:) = tmphdxf(n,:)
    tmpqc(nn) = tmpqc(n)
  END DO
  nobs = nn
  WRITE(6,'(I10,A,I3.3)') nobs,' OBSERVATIONS TO BE ASSIMILATED IN MYRANK ',myrank
!
! SORT
!
  ALLOCATE( tmp2elm(nobs) )
  ALLOCATE( tmp2lon(nobs) )
  ALLOCATE( tmp2lat(nobs) )
  ALLOCATE( tmp2lev(nobs) )
  ALLOCATE( tmp2dat(nobs) )
  ALLOCATE( tmp2err(nobs) )
!  ALLOCATE( tmp2k(nobs) )
  ALLOCATE( tmp2dep(nobs) )
  ALLOCATE( tmp2hdxf(nobs,nbv) )
  ALLOCATE( obselm(nobs) )
  ALLOCATE( obslon(nobs) )
  ALLOCATE( obslat(nobs) )
  ALLOCATE( obslev(nobs) )
  ALLOCATE( obsdat(nobs) )
  ALLOCATE( obserr(nobs) )
!  ALLOCATE( obsk(nobs) )
  ALLOCATE( obsdep(nobs) )
  ALLOCATE( obshdxf(nobs,nbv) )
  nobsgrd = 0
  nj = 0
!$OMP PARALLEL PRIVATE(i,j,n,nn)
!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    DO n=1,nobs
      IF(tmplat(n) < lat(j) .OR. lat(j+1) <= tmplat(n)) CYCLE
      nj(j) = nj(j) + 1
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    njs(j) = SUM(nj(0:j-1))
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    nn = 0
    DO n=1,nobs
      IF(tmplat(n) < lat(j) .OR. lat(j+1) <= tmplat(n)) CYCLE
      nn = nn + 1
      tmp2elm(njs(j)+nn) = tmpelm(n)
      tmp2lon(njs(j)+nn) = tmplon(n)
      tmp2lat(njs(j)+nn) = tmplat(n)
      tmp2lev(njs(j)+nn) = tmplev(n)
      tmp2dat(njs(j)+nn) = tmpdat(n)
      tmp2err(njs(j)+nn) = tmperr(n)
!      tmp2k(njs(j)+nn) = tmpk(n)
      tmp2dep(njs(j)+nn) = tmpdep(n)
      tmp2hdxf(njs(j)+nn,:) = tmphdxf(n,:)
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC)
  DO j=1,nlat-1
    IF(nj(j) == 0) THEN
      nobsgrd(:,j) = njs(j)
      CYCLE
    END IF
    nn = 0
    DO i=1,nlon
      DO n=njs(j)+1,njs(j)+nj(j)
        IF(i < nlon) THEN
          IF(tmp2lon(n) < lon(i) .OR. lon(i+1) <= tmp2lon(n)) CYCLE
        ELSE
          IF(tmp2lon(n) < lon(nlon) .OR. 360.0d0 <= tmp2lon(n)) CYCLE
        END IF
        nn = nn + 1
        obselm(njs(j)+nn) = tmp2elm(n)
        obslon(njs(j)+nn) = tmp2lon(n)
        obslat(njs(j)+nn) = tmp2lat(n)
        obslev(njs(j)+nn) = tmp2lev(n)
        obsdat(njs(j)+nn) = tmp2dat(n)
        obserr(njs(j)+nn) = tmp2err(n)
!        obsk(njs(j)+nn) = tmp2k(n)
        obsdep(njs(j)+nn) = tmp2dep(n)
        obshdxf(njs(j)+nn,:) = tmp2hdxf(n,:)
      END DO
      nobsgrd(i,j) = njs(j) + nn
    END DO
    IF(nn /= nj(j)) THEN
!$OMP CRITICAL
      WRITE(6,'(A,2I)') 'OBS DATA SORT ERROR: ',nn,nj(j)
      WRITE(6,'(F6.2,A,F6.2)') lat(j),'< LAT <',lat(j+1)
      WRITE(6,'(F6.2,A,F6.2)') MINVAL(tmp2lat(njs(j)+1:njs(j)+nj(j))),'< OBSLAT <',MAXVAL(tmp2lat(njs(j)+1:njs(j)+nj(j)))
!$OMP END CRITICAL
    END IF
  END DO
!$OMP END DO
!$OMP END PARALLEL
  DEALLOCATE( tmp2elm )
  DEALLOCATE( tmp2lon )
  DEALLOCATE( tmp2lat )
  DEALLOCATE( tmp2lev )
  DEALLOCATE( tmp2dat )
  DEALLOCATE( tmp2err )
!  DEALLOCATE( tmp2k )
  DEALLOCATE( tmp2dep )
  DEALLOCATE( tmp2hdxf )
  DEALLOCATE( tmpelm )
  DEALLOCATE( tmplon )
  DEALLOCATE( tmplat )
  DEALLOCATE( tmplev )
  DEALLOCATE( tmpdat )
  DEALLOCATE( tmperr )
  DEALLOCATE( tmpk )
  DEALLOCATE( tmpdep )
  DEALLOCATE( tmphdxf )
  DEALLOCATE( tmpqc )

  RETURN
END SUBROUTINE set_letkf_obs
!-----------------------------------------------------------------------
! Monitor departure from gues/anal mean
!-----------------------------------------------------------------------
SUBROUTINE monit_mean(file)
  IMPLICIT NONE
  CHARACTER(4),INTENT(IN) :: file
  REAL(r_size) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size) :: v2d(nlon,nlat,nv2d)
  REAL(r_size) :: p_full(nlon,nlat,nlev)
  REAL(r_size) :: elem
  REAL(r_size) :: bias_u,bias_v,bias_t,bias_ps,bias_q,bias_rh
  REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_ps,rmse_q,rmse_rh
  REAL(r_size) :: hdxf,dep,ri,rj,rk
  INTEGER :: n,iu,iv,it,iq,ips,irh
  CHARACTER(11) :: filename='filexxx.grd'

  rmse_u  = 0.0d0
  rmse_v  = 0.0d0
  rmse_t  = 0.0d0
  rmse_q  = 0.0d0
  rmse_ps = 0.0d0
  rmse_rh = 0.0d0
  bias_u = 0.0d0
  bias_v = 0.0d0
  bias_t = 0.0d0
  bias_q = 0.0d0
  bias_ps = 0.0d0
  bias_rh = 0.0d0
  iu  = 0
  iv  = 0
  it  = 0
  iq  = 0
  ips = 0
  irh = 0

  WRITE(filename(1:7),'(A4,A3)') file,'_me'
  CALL read_grd(filename,v3d,v2d)
  CALL calc_pfull(nlon,nlat,v2d(:,:,iv2d_ps),p_full)

  DO n=1,nobs
    CALL phys2ijk(p_full,obselm(n),obslon(n),obslat(n),obslev(n),ri,rj,rk)
    IF(CEILING(rk) > nlev) CYCLE
    IF(CEILING(rk) < 2 .AND. NINT(obselm(n)) /= id_ps_obs) THEN
      IF(NINT(obselm(n)) == id_u_obs .OR. NINT(obselm(n)) == id_v_obs) THEN
        rk = 1.00001d0
      ELSE
        CYCLE
      END IF
    END IF
    CALL Trans_XtoY(obselm(n),ri,rj,rk,v3d,v2d,p_full,hdxf)
    dep = obsdat(n) - hdxf
    SELECT CASE(NINT(obselm(n)))
    CASE(id_u_obs)
      rmse_u = rmse_u + dep**2
      bias_u = bias_u + dep
      iu = iu + 1
    CASE(id_v_obs)
      rmse_v = rmse_v + dep**2
      bias_v = bias_v + dep
      iv = iv + 1
    CASE(id_t_obs)
      rmse_t = rmse_t + dep**2
      bias_t = bias_t + dep
      it = it + 1
    CASE(id_q_obs)
      rmse_q = rmse_q + dep**2
      bias_q = bias_q + dep
      iq = iq + 1
    CASE(id_ps_obs)
      rmse_ps = rmse_ps + dep**2
      bias_ps = bias_ps + dep
      ips = ips + 1
    CASE(id_rh_obs)
      rmse_rh = rmse_rh + dep**2
      bias_rh = bias_rh + dep
      irh = irh + 1
    END SELECT
  END DO

  IF(iu == 0) THEN
    rmse_u = undef
    bias_u = undef
  ELSE
    rmse_u = SQRT(rmse_u / REAL(iu,r_size))
    bias_u = bias_u / REAL(iu,r_size)
  END IF
  IF(iv == 0) THEN
    rmse_v = undef
    bias_v = undef
  ELSE
    rmse_v = SQRT(rmse_v / REAL(iv,r_size))
    bias_v = bias_v / REAL(iv,r_size)
  END IF
  IF(it == 0) THEN
    rmse_t = undef
    bias_t = undef
  ELSE
    rmse_t = SQRT(rmse_t / REAL(it,r_size))
    bias_t = bias_t / REAL(it,r_size)
  END IF
  IF(iq == 0) THEN
    rmse_q = undef
    bias_q = undef
  ELSE
    rmse_q = SQRT(rmse_q / REAL(iq,r_size))
    bias_q = bias_q / REAL(iq,r_size)
  END IF
  IF(ips == 0) THEN
    rmse_ps = undef
    bias_ps = undef
  ELSE
    rmse_ps = SQRT(rmse_ps / REAL(ips,r_size))
    bias_ps = bias_ps / REAL(ips,r_size)
  END IF
  IF(irh == 0) THEN
    rmse_rh = undef
    bias_rh = undef
  ELSE
    rmse_rh = SQRT(rmse_rh / REAL(irh,r_size))
    bias_rh = bias_rh / REAL(irh,r_size)
  END IF

  WRITE(6,'(3A)') '== PARTIAL OBSERVATIONAL DEPARTURE (',file,') =============================='
  WRITE(6,'(6A12)') 'U','V','T','Q','PS','RH'
  WRITE(6,'(6ES12.3)') bias_u,bias_v,bias_t,bias_q,bias_ps,bias_rh
  WRITE(6,'(6ES12.3)') rmse_u,rmse_v,rmse_t,rmse_q,rmse_ps,rmse_rh
  WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS =============================================='
  WRITE(6,'(6A12)') 'U','V','T','Q','PS','RH'
  WRITE(6,'(6I12)') iu,iv,it,iq,ips,irh
  WRITE(6,'(A)') '========================================================================'

  RETURN
END SUBROUTINE monit_mean
!------------------------------------------------------------------------
! search a latitude with the longest lon-grids for the localization mask
!------------------------------------------------------------------------
SUBROUTINE search_longestlat(dist_zero,lat_cnt,dlat_zero,dlon_local)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN)  :: dist_zero
  REAL(r_size),INTENT(IN)  :: lat_cnt
  REAL(r_size),INTENT(IN)  :: dlat_zero
  REAL(r_size),INTENT(OUT) :: dlon_local
  REAL(r_size) :: lat_dist,lat_c
  REAL(r_size) :: cosd
  INTEGER :: j,j0,j1

  IF(ABS(lat_cnt)+dlat_zero >= 90.0d0) THEN
    dlon_local = 360.0d0
!   dlon_local = 180.0d0
    RETURN
  ENDIF

!DEBUG by kotsuki!   DO j=1,nlat
!DEBUG by kotsuki!    IF( lat(j)>=(lat_cnt-dlat_zero) .AND. lat(j-1)<(lat_cnt-dlat_zero) ) THEN
!DEBUG by kotsuki!      j0 = j
!DEBUG by kotsuki!    ENDIF
!DEBUG by kotsuki!    IF( lat(j)>=(lat_cnt+dlat_zero) ) THEN
!DEBUG by kotsuki!      j1 = j-1
!DEBUG by kotsuki!      EXIT
!DEBUG by kotsuki!    ENDIF
!DEBUG by kotsuki!  END DO

  lat_c = lat_cnt * pi/180.0d0
  dlon_local = 0.0d0
!DEBUG by kotsuki!  DO j=j0,j1
  DO j=1,nlat
    lat_dist = lat(j) * pi/180.0d0
    cosd = (COS(dist_zero/re)-SIN(lat_dist)*SIN(lat_c)) &
   &       / (COS(lat_dist)*COS(lat_c))
    cosd = MIN(1.0d0,cosd)
    cosd = MAX(-1.0d0,cosd)
    dlon_local = MAX(ACOS(cosd)/pi*180.0d0,dlon_local)
  END DO

  RETURN
END SUBROUTINE search_longestlat
END MODULE letkf_obs
