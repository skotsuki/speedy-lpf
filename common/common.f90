MODULE common
!=======================================================================
!
! [PURPOSE:] General constants and procedures
!
! [ATTENTION:] This module calls 'SFMT.f90'
!
! [HISTORY:]
!   07/20/2004 Takemasa MIYOSHI  created
!   01/23/2009 Takemasa MIYOSHI  modified for SFMT
!
!=======================================================================
  IMPLICIT NONE
  PUBLIC
  INTEGER,SAVE :: nbv=02
!-----------------------------------------------------------------------
! Variable size definitions
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
!-----------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------
  REAL(r_size),PARAMETER :: pi=3.1415926535d0
  REAL(r_size),PARAMETER :: gg=9.81d0
  REAL(r_size),PARAMETER :: rd=287.05d0       ! gas constant air (J/kg/K)      GYL
  REAL(r_size),PARAMETER :: rv=461.50d0       ! gas constant H2O (J/kg/K)      GYL
  REAL(r_size),PARAMETER :: cp=1005.7d0       ! spec heat air [p] (J/kg/K)
  REAL(r_size),PARAMETER :: hvap=2.5d6        ! heat of vaporization (J/kg)
  REAL(r_size),PARAMETER :: fvirt=rv/rd-1.0d0 ! parameter for T/Tv conversion  GYL
  REAL(r_size),PARAMETER :: re=6371.3d3
  REAL(r_size),PARAMETER :: r_omega=7.292d-5
  REAL(r_size),PARAMETER :: t0c=273.15d0
  REAL(r_size),PARAMETER :: undef=-9.99d33

CONTAINS
!-----------------------------------------------------------------------
! Mean
!-----------------------------------------------------------------------
SUBROUTINE com_mean(ndim,var,amean)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var(ndim)
  REAL(r_size),INTENT(OUT) :: amean

  INTEGER :: i

  amean = 0.0d0
  DO i=1,ndim
    amean = amean + var(i)
  END DO
  amean = amean / REAL(ndim,r_size)

  RETURN
END SUBROUTINE com_mean
!-----------------------------------------------------------------------
! Standard deviation
!-----------------------------------------------------------------------
SUBROUTINE com_stdev(ndim,var,aout)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var(ndim)
  REAL(r_size),INTENT(OUT) :: aout

  REAL(r_size) :: amean
  REAL(r_size) :: dev(ndim)

  CALL com_mean(ndim,var,amean)

  dev(:) = var(:) - amean

  aout = SQRT( SUM(dev*dev) / REAL(ndim-1,r_size) )

  RETURN
END SUBROUTINE com_stdev
!-----------------------------------------------------------------------
! Covariance
!-----------------------------------------------------------------------
SUBROUTINE com_covar(ndim,var1,var2,cov)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var1(ndim)
  REAL(r_size),INTENT(IN) :: var2(ndim)
  REAL(r_size),INTENT(OUT) :: cov

  REAL(r_size) :: amean1,amean2
  REAL(r_size) :: dev1(ndim),dev2(ndim)

  CALL com_mean(ndim,var1,amean1)
  CALL com_mean(ndim,var2,amean2)

  dev1(:) = var1(:) - amean1
  dev2(:) = var2(:) - amean2

  cov = SUM( dev1*dev2 ) / REAL(ndim-1,r_size)

  RETURN
END SUBROUTINE com_covar
!-----------------------------------------------------------------------
! Correlation
!-----------------------------------------------------------------------
SUBROUTINE com_correl(ndim,var1,var2,cor)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var1(ndim)
  REAL(r_size),INTENT(IN) :: var2(ndim)
  REAL(r_size),INTENT(OUT) :: cor

  REAL(r_size) :: cov,stdev1,stdev2

  CALL com_stdev(ndim,var1,stdev1)
  CALL com_stdev(ndim,var2,stdev2)
  CALL com_covar(ndim,var1,var2,cov)

  cor = cov/stdev1/stdev2

  RETURN
END SUBROUTINE com_correl
!-----------------------------------------------------------------------
! Anomaly Correlation
!-----------------------------------------------------------------------
SUBROUTINE com_anomcorrel(ndim,var1,var2,varmean,cor)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var1(ndim)
  REAL(r_size),INTENT(IN) :: var2(ndim)
  REAL(r_size),INTENT(IN) :: varmean(ndim)
  REAL(r_size),INTENT(OUT) :: cor

  REAL(r_size) :: dev1(ndim),dev2(ndim)

  dev1 = var1 - varmean
  dev2 = var2 - varmean

  cor = SUM( dev1*dev2 ) / SQRT( SUM(dev1*dev1) * SUM(dev2*dev2) )

  RETURN
END SUBROUTINE com_anomcorrel
!-----------------------------------------------------------------------
! L2 Norm
!-----------------------------------------------------------------------
SUBROUTINE com_l2norm(ndim,var,anorm)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var(ndim)
  REAL(r_size),INTENT(OUT) :: anorm

  anorm = SQRT( SUM(var*var) )

  RETURN
END SUBROUTINE com_l2norm
!-----------------------------------------------------------------------
! RMS (root mean square)
!-----------------------------------------------------------------------
SUBROUTINE com_rms(ndim,var,rmsv)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var(ndim)
  REAL(r_size),INTENT(OUT) :: rmsv

  rmsv = SQRT( SUM(var*var) / REAL(ndim,r_size) )

  RETURN
END SUBROUTINE com_rms
!-----------------------------------------------------------------------
! Lanczos Filter (Low-pass) with cyclic boundary
!-----------------------------------------------------------------------
SUBROUTINE com_filter_lanczos(ndim,fc,var)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: fc    ! critical frequency in [0,pi]
  REAL(r_size),INTENT(INOUT) :: var(ndim)

  INTEGER,PARAMETER :: lresol=10

  REAL(r_size) :: weight(-lresol:lresol)
  REAL(r_size) :: varwk(1-lresol:ndim+lresol)
  REAL(r_size) :: rl,rlresol
  INTEGER :: i,l
!
! Weight
!
  rlresol = REAL(lresol,r_size)
  DO l=-lresol,-1
    rl = REAL(l,r_size)
    weight(l) = SIN(fc*rl) * SIN(pi*rl/rlresol) &
      & * rlresol / pi / rl / pi / rl
  END DO
  DO l=1,lresol
    rl = REAL(l,r_size)
    weight(l) = SIN(fc*rl) * SIN(pi*rl/rlresol) &
      & * rlresol / pi / rl / pi / rl
  END DO
  weight(0) = fc / pi
!
! Cyclic boundary
!
  DO i=0,1-lresol,-1
    varwk(i) = var(ndim+i)
  END DO
  DO i=ndim+1,ndim+lresol
    varwk(i) = var(i-ndim)
  END DO
  varwk(1:ndim) = var(1:ndim)
!
! Filter
!
  var = 0.0d0
  DO i=1,ndim
    DO l=-lresol,lresol
      var(i) = var(i) + weight(l) * varwk(i+l)
    END DO
  END DO

  RETURN
END SUBROUTINE com_filter_lanczos
!-----------------------------------------------------------------------
! RAND (random number with uniform distribution)
!-----------------------------------------------------------------------
SUBROUTINE com_rand(ndim,var)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(OUT) :: var(1:ndim)
  REAL(r_dble) :: genrand_res53
  INTEGER :: idate(8)
  INTEGER :: i,iseed
  LOGICAL,SAVE :: first=.true.

  IF (first) THEN
    CALL DATE_AND_TIME(VALUES=idate)
    iseed = idate(8) + idate(7)*1000
    !iseed = 1 ! for debug, SK 20190425 for LPF to reproduce the same experiments
    CALL init_gen_rand(iseed)
    first=.false.
  END IF

  DO i=1,ndim
    var(i) = genrand_res53()
  END DO

  RETURN
END SUBROUTINE com_rand
!-----------------------------------------------------------------------
SUBROUTINE com_rand_seed(ndim,iseed,var) ! from KK 20200426
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim, iseed
  REAL(r_size),INTENT(OUT) :: var(1:ndim)
  REAL(r_dble) :: genrand_res53
  INTEGER :: idate(8)
  INTEGER :: i, jseed
  LOGICAL,SAVE :: first=.true.

  IF (first) THEN
    !!!print *, first, iseed
    IF( iseed==0 ) THEN
      CALL DATE_AND_TIME(VALUES=idate)
      jseed = idate(8) + idate(7)*1000
    ELSE
      jseed = iseed
    ENDIF
    CALL init_gen_rand(jseed)
    first=.false.
  END IF

  DO i=1,ndim
    var(i) = genrand_res53()
  END DO

  RETURN
END SUBROUTINE com_rand_seed
!-----------------------------------------------------------------------
! RANDN (random number with normal distribution)
!-----------------------------------------------------------------------
SUBROUTINE com_randn(ndim,var)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(OUT) :: var(1:ndim)
  REAL(r_size) :: rnd(2)
  REAL(r_dble) :: genrand_res53
  INTEGER :: idate(8)
  INTEGER :: i,iseed
  LOGICAL,SAVE :: first=.true.

  IF (first) THEN
    CALL DATE_AND_TIME(VALUES=idate)
    iseed = idate(8) + idate(7)*1000
    CALL init_gen_rand(iseed)
    first=.false.
  END IF

  IF( MOD(ndim,2)==0 ) THEN
    DO i=1,ndim/2
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
    END DO
  ELSE
    DO i=1,(ndim-1)/2
      rnd(1) = genrand_res53()
      rnd(2) = genrand_res53()
      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
    END DO
    rnd(1) = genrand_res53()
    rnd(2) = genrand_res53()
    var(ndim) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
  END IF

  RETURN
END SUBROUTINE com_randn
!NO,SOFAR!!-----------------------------------------------------------------------
!NO,SOFAR!SUBROUTINE com_randn_seed(ndim,iseed,var) ! from KK 20200426
!NO,SOFAR!  IMPLICIT NONE
!NO,SOFAR!
!NO,SOFAR!  INTEGER,INTENT(IN) :: ndim, iseed
!NO,SOFAR!  REAL(r_size),INTENT(OUT) :: var(1:ndim)
!NO,SOFAR!  REAL(r_size) :: rnd(2)
!NO,SOFAR!  REAL(r_dble) :: genrand_res53
!NO,SOFAR!  INTEGER :: i
!NO,SOFAR!  LOGICAL,SAVE :: first=.true.
!NO,SOFAR!
!NO,SOFAR!  IF (first) THEN
!NO,SOFAR!    CALL init_gen_rand(iseed)
!NO,SOFAR!    first=.true.
!NO,SOFAR!  END IF
!NO,SOFAR!
!NO,SOFAR!  IF( MOD(ndim,2)==0 ) THEN
!NO,SOFAR!    DO i=1,ndim/2
!NO,SOFAR!      rnd(1) = genrand_res53()
!NO,SOFAR!      rnd(2) = genrand_res53()
!NO,SOFAR!      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
!NO,SOFAR!      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
!NO,SOFAR!    END DO
!NO,SOFAR!  ELSE
!NO,SOFAR!    DO i=1,(ndim-1)/2
!NO,SOFAR!      rnd(1) = genrand_res53()
!NO,SOFAR!      rnd(2) = genrand_res53()
!NO,SOFAR!      var(i*2-1) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
!NO,SOFAR!      var(i*2) = sqrt( -2.0d0 * log( rnd(1) ) ) * cos( 2.0d0*pi*rnd(2) )
!NO,SOFAR!    END DO
!NO,SOFAR!    rnd(1) = genrand_res53()
!NO,SOFAR!    rnd(2) = genrand_res53()
!NO,SOFAR!    var(ndim) = sqrt( -2.0d0 * log( rnd(1) ) ) * sin( 2.0d0*pi*rnd(2) )
!NO,SOFAR!  END IF
!NO,SOFAR!
!NO,SOFAR!  RETURN
!NO,SOFAR!END SUBROUTINE com_randn_seed
!-----------------------------------------------------------------------
! TIMEINC
!-----------------------------------------------------------------------
SUBROUTINE com_timeinc_hr(iy,im,id,ih,incr)
  IMPLICIT NONE

  INTEGER,INTENT(INOUT) :: iy
  INTEGER,INTENT(INOUT) :: im
  INTEGER,INTENT(INOUT) :: id
  INTEGER,INTENT(INOUT) :: ih
  INTEGER,INTENT(IN) :: incr

  ih = ih + incr
  IF(ih>23) THEN
    ih = ih - 24
    id = id + 1
    IF(id==29.AND.im==2.AND.mod(iy,4)/=0) THEN
      id = 1
      im = 3
    ELSE IF(id==30.AND.im==2.AND.mod(iy,4)==0) THEN
      id = 1
      im = 3
    ELSE IF(id==31.AND.(im==4.OR.im==6.OR.im==9.OR.im==11)) THEN
      id = 1
      im = im + 1
    ELSE IF(id==32.AND.(im==1.OR.im==3.OR.im==5.OR.im==7.OR.im==8.OR.im==10)) THEN
      id = 1
      im = im + 1
    ELSE IF(id==32.AND.im==12) THEN
      id = 1
      im = 1
      iy = iy + 1
    END IF
  END IF

  RETURN
END SUBROUTINE com_timeinc_hr
!-----------------------------------------------------------------------
! TIMECONVERSION
!-----------------------------------------------------------------------
SUBROUTINE com_time2ymdh(itime,iy,im,id,ih)
  IMPLICIT NONE
  INTEGER(8),INTENT(IN) :: itime
  INTEGER,INTENT(OUT) :: iy
  INTEGER,INTENT(OUT) :: im
  INTEGER,INTENT(OUT) :: id
  INTEGER,INTENT(OUT) :: ih

  iy = INT(  itime / 1000000 )
  im = INT( (itime-iy*1000000) / 10000 )
  id = INT( (itime-iy*1000000-im*10000) / 100 )
  ih = INT(  itime-iy*1000000-im*10000-id*100 )

  RETURN
END SUBROUTINE com_time2ymdh

SUBROUTINE com_ymdh2time(iy,im,id,ih,itime)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: iy
  INTEGER,INTENT(IN) :: im
  INTEGER,INTENT(IN) :: id
  INTEGER,INTENT(IN) :: ih
  INTEGER(8),INTENT(OUT) :: itime

  itime=iy*1000000+im*10000+id*100+ih

  RETURN
END SUBROUTINE com_ymdh2time
!-----------------------------------------------------------------------
! DISTANCE BETWEEN TWO POINTS (LONa,LATa)-(LONb,LATb)
!-----------------------------------------------------------------------
SUBROUTINE com_distll(ndim,alon,alat,blon,blat,dist)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: alon(ndim)
  REAL(r_size),INTENT(IN) :: alat(ndim)
  REAL(r_size),INTENT(IN) :: blon(ndim)
  REAL(r_size),INTENT(IN) :: blat(ndim)
  REAL(r_size),INTENT(OUT) :: dist(ndim)
  REAL(r_size),PARAMETER :: r180=1.0d0/180.0d0
  REAL(r_size) :: lon1,lon2,lat1,lat2
  REAL(r_size) :: cosd(ndim)
  INTEGER :: i

  DO i=1,ndim
    lon1 = alon(i) * pi * r180
    lon2 = blon(i) * pi * r180
    lat1 = alat(i) * pi * r180
    lat2 = blat(i) * pi * r180

    cosd(i) = SIN(lat1)*SIN(lat2) + COS(lat1)*COS(lat2)*COS(lon2-lon1)
    cosd(i) = MIN( 1.d0,cosd(i))
    cosd(i) = MAX(-1.d0,cosd(i))

    dist(i) = ACOS( cosd(i) ) * re
  END DO

  RETURN
END SUBROUTINE com_distll
!-----------------------------------------------------------------------
! DISTANCE BETWEEN TWO POINTS (LONa,LATa)-(LONb,LATb)
!-----------------------------------------------------------------------
SUBROUTINE com_distll_1(alon,alat,blon,blat,dist)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: alon
  REAL(r_size),INTENT(IN) :: alat
  REAL(r_size),INTENT(IN) :: blon
  REAL(r_size),INTENT(IN) :: blat
  REAL(r_size),INTENT(OUT) :: dist
  REAL(r_size),PARAMETER :: r180=1.0d0/180.0d0
  REAL(r_size) :: lon1,lon2,lat1,lat2
  REAL(r_size) :: cosd

  lon1 = alon * pi * r180
  lon2 = blon * pi * r180
  lat1 = alat * pi * r180
  lat2 = blat * pi * r180

  cosd = SIN(lat1)*SIN(lat2) + COS(lat1)*COS(lat2)*COS(lon2-lon1)
  cosd = MIN( 1.d0,cosd)
  cosd = MAX(-1.d0,cosd)

  dist = ACOS( cosd ) * re

  RETURN
END SUBROUTINE com_distll_1

FUNCTION com_distll_1f(alon,alat,blon,blat) RESULT(dist)
  IMPLICIT NONE
  REAL(r_size) :: alon
  REAL(r_size) :: alat
  REAL(r_size) :: blon
  REAL(r_size) :: blat
  REAL(r_size) :: dist
  REAL(r_size),PARAMETER :: r180=1.0d0/180.0d0
  REAL(r_size) :: lon1,lon2,lat1,lat2
  REAL(r_size) :: cosd

  lon1 = alon * pi * r180
  lon2 = blon * pi * r180
  lat1 = alat * pi * r180
  lat2 = blat * pi * r180

  cosd = SIN(lat1)*SIN(lat2) + COS(lat1)*COS(lat2)*COS(lon2-lon1)
  cosd = MIN( 1.d0,cosd)
  cosd = MAX(-1.d0,cosd)

  dist = ACOS( cosd ) * re

  RETURN
END FUNCTION com_distll_1f
!-----------------------------------------------------------------------
! Cubic spline interpolation
!   [Reference:] Akima, H., 1970: J. ACM, 17, 589-602.
!-----------------------------------------------------------------------
SUBROUTINE com_interp_spline(ndim,x,y,n,x5,y5)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ndim         ! number of grid points
  REAL(r_size),INTENT(IN) :: x(ndim) ! coordinate
  REAL(r_size),INTENT(IN) :: y(ndim) ! variable
  INTEGER,INTENT(IN) :: n            ! number of targets
  REAL(r_size),INTENT(IN) :: x5(n)   ! target coordinates
  REAL(r_size),INTENT(OUT) :: y5(n)  ! target values
  INTEGER :: i,j,m
  REAL(r_size) :: dydx(5),ddydx(4),t(2),dx21,dx
  REAL(r_size) :: wk

  TGT: DO j=1,n
    DO i=1,ndim
      IF(x5(j) == x(i)) THEN
        y5(j) = y(i)
        CYCLE TGT
      END IF
      IF(x5(j) < x(i)) EXIT
    END DO
!       i-3   i-2   i-1    i    i+1   i+2
!     ---+-----+-----+---*-+-----+-----+---
!dydx       1     2     3     4     5
!ddydx         1     2     3     4
!t                   1     2
    IF(i==2) THEN
      DO m=3,5
        dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
      END DO
      dydx(2) = 2.0d0*dydx(3) - dydx(4)
      dydx(1) = 2.0d0*dydx(2) - dydx(3)
    ELSE IF(i==3) THEN
      DO m=2,5
        dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
      END DO
      dydx(1) = 2.0d0*dydx(2) - dydx(3)
    ELSE IF(i==ndim) THEN
      DO m=1,3
        dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
      END DO
      dydx(4) = 2.0d0*dydx(3) - dydx(2)
      dydx(5) = 2.0d0*dydx(4) - dydx(3)
    ELSE IF(i==ndim-1) THEN
      DO m=1,4
        dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
      END DO
      dydx(5) = 2.0d0*dydx(4) - dydx(3)
    ELSE
      DO m=1,5
        dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
      END DO
    END IF
    DO m=1,4
      ddydx(m) = ABS(dydx(m+1) - dydx(m))
    END DO
    DO m=1,2
      wk = ddydx(m+2) + ddydx(m)
      IF(wk == 0) THEN
        t(m) = 0.0d0
      ELSE
        t(m) = (ddydx(m+2)*dydx(m+1)+ddydx(m)*dydx(m+2))/wk
      END IF
    END DO
    dx21 = x(i)-x(i-1)
    dx = x5(j) - x(i-1)
    y5(j) = y(i-1) &
        & + dx*t(1) &
        & + dx*dx*(3.0d0*dydx(3)-2.0d0*t(1)-t(2))/dx21 &
        & + dx*dx*dx*(t(1)+t(2)-2.0d0*dydx(3))/dx21/dx21
  END DO TGT

  RETURN
END SUBROUTINE com_interp_spline
!-----------------------------------------------------------------------
! (LON,LAT) --> (i,j) conversion
!   [ORIGINAL AUTHOR:] Masaru Kunii
!-----------------------------------------------------------------------
SUBROUTINE com_pos2ij(msw,nx,ny,flon,flat,num_obs,olon,olat,oi,oj)
  IMPLICIT NONE
  ! --- inout variables
  INTEGER,INTENT(IN) :: msw   !MODE SWITCH: 1: fast, 2: accurate
  INTEGER,INTENT(IN) :: nx,ny !number of grid points
  REAL(r_size),INTENT(IN) :: flon(nx,ny),flat(nx,ny) !(lon,lat) at (i,j)
  INTEGER,INTENT(IN) :: num_obs !repetition number of conversion
  REAL(r_size),INTENT(IN) :: olon(num_obs),olat(num_obs) !target (lon,lat)
  REAL(r_size),INTENT(OUT) :: oi(num_obs),oj(num_obs) !target (i,j)
  ! --- local work variables
  LOGICAL,PARAMETER :: detailout = .FALSE.
  INTEGER,PARAMETER :: num_grid_ave = 4  ! fix
  INTEGER :: inum,ix,jy,ip,wk_maxp
  INTEGER :: iorder_we,iorder_sn
  INTEGER :: nxp,nyp
  REAL(r_size),PARAMETER :: miss = -32768 
  REAL(r_size),PARAMETER :: max_dist = 2.0e+6
  REAL(r_size) :: rlat_max, rlat_min, rlon_max, rlon_min   
  REAL(r_size) :: dist(num_grid_ave) 
  REAL(r_size) :: dist_min_x(num_obs, num_grid_ave)
  REAL(r_size) :: dist_min_y(num_obs, num_grid_ave) 
  REAL(r_size) :: wk_dist, sum_dist
  REAL(r_size) :: ratio(num_grid_ave)
  IF(detailout) THEN
    WRITE(6,'(A)') '====================================================='
    WRITE(6,'(A)') '      Detailed output of SUBROUTINE com_pos2ij       '
    WRITE(6,'(A)') '====================================================='    
  END IF
  ! ================================================================
  !   Check the Order of flon, flat
  ! ================================================================   
  iorder_we = 1
  iorder_sn = 1
  IF(flon(1,1) > flon(2,1)) THEN
    iorder_we = -1
  END IF
  IF(flat(1,1) > flat(1,2)) THEN
    iorder_sn = -1
  END IF
  IF(detailout) THEN  
    WRITE(6,'(3X,A,I5)') 'Obs Order (WE) :',iorder_we 
    WRITE(6,'(3X,A,I5)') 'Obs Order (SN) :',iorder_sn 
  END IF
  ! ================================================================
  !  FAST MODE
  ! ================================================================   
  IF(msw == 1) THEN
    ! ==============================================================
    !   Surrounding 4 Grid Points Interpolation
    ! ==============================================================   
    Obs_Loop_1 : DO inum=1,num_obs 
      IF(detailout) WRITE(6,'(A,I5,2F15.5)') '*** START OBS ',inum,olon(inum),olat(inum) 
      ! ------------------------------------------------------------
      !    Search Basic Point
      ! ------------------------------------------------------------ 
      nxp = miss
      nyp = miss
      DO jy=1,ny-1
        DO ix=1,nx-1
          rlon_max = MAXVAL(flon(ix:ix+1, jy:jy+1))
          rlon_min = MINVAL(flon(ix:ix+1, jy:jy+1))
          rlat_max = MAXVAL(flat(ix:ix+1, jy:jy+1))
          rlat_min = MINVAL(flat(ix:ix+1, jy:jy+1))
          IF(rlon_min <= olon(inum) .AND. rlon_max >= olon(inum) .AND. &
           & rlat_min <= olat(inum) .AND. rlat_max >= olat(inum)) THEN
            nxp = ix
            nyp = jy
            EXIT
          END IF
        END DO
      END DO
      IF(detailout) WRITE(6,'(3X,A,2I7)') 'nxp, nyp =',nxp,nyp
      IF(nxp == miss .OR. nyp == miss) THEN
        WRITE(6,'(A)') '!!WARNING(com_pos2ij): obs position cannot be detected'
        oi(inum) = miss
        oj(inum) = miss
        CYCLE Obs_Loop_1
      END IF
      ! ------------------------------------------------------------
      !    Interpolation
      ! ------------------------------------------------------------    
      CALL com_distll_1(flon(nxp  ,nyp  ),flat(nxp  ,nyp  ),&
                      & olon(inum),olat(inum),dist(1))
      CALL com_distll_1(flon(nxp+1,nyp  ),flat(nxp+1,nyp  ),&
                      & olon(inum),olat(inum),dist(2))
      CALL com_distll_1(flon(nxp  ,nyp+1),flat(nxp  ,nyp+1),&
                      & olon(inum),olat(inum),dist(3))      
      CALL com_distll_1(flon(nxp+1,nyp+1),flat(nxp+1,nyp+1),&
                      & olon(inum),olat(inum),dist(4))      
      dist(1:4) = dist(1:4) * 1.D-3  
      IF(detailout) WRITE(6,'(3X,A,4F15.5)') 'distance :',dist(1:4) 
      sum_dist = dist(1) * dist(1) * dist(2) * dist(2) * dist(3) * dist(3) &
             & + dist(2) * dist(2) * dist(3) * dist(3) * dist(4) * dist(4) &
             & + dist(3) * dist(3) * dist(4) * dist(4) * dist(1) * dist(1) &
             & + dist(4) * dist(4) * dist(1) * dist(1) * dist(2) * dist(2)
      ratio(1) = (dist(2)*dist(2)*dist(3)*dist(3)*dist(4)*dist(4))/sum_dist
      ratio(2) = (dist(3)*dist(3)*dist(4)*dist(4)*dist(1)*dist(1))/sum_dist
      ratio(3) = (dist(4)*dist(4)*dist(1)*dist(1)*dist(2)*dist(2))/sum_dist
      ratio(4) = (dist(1)*dist(1)*dist(2)*dist(2)*dist(3)*dist(3))/sum_dist
      IF(detailout) WRITE(6,'(3X,A,5F15.5)') 'ratio    :',ratio(1:4),SUM(ratio(1:4))
      oi(inum) = ratio(1) *  nxp    + ratio(2) * (nxp+1) &
             & + ratio(3) *  nxp    + ratio(4) * (nxp+1)
      oj(inum) = ratio(1) *  nyp    + ratio(2) *  nyp    &
             & + ratio(3) * (nyp+1) + ratio(4) * (nyp+1)    
      IF(detailout) WRITE(6,'(3X,A,2F15.5)') 'position :',oi(inum), oj(inum)
 
    END DO Obs_Loop_1
  ! ================================================================
  !  ACCURATE MODE
  ! ================================================================   
  ELSE IF(msw == 2) THEN
    ! ================================================================
    !   Nearest 4 Grid Points Interpolation
    ! ================================================================   
    Obs_Loop_2 : DO inum=1,num_obs
      IF(detailout) WRITE(6,'(A,I5,2F15.5)') '*** START OBS ',inum,olon(inum),olat(inum) 
      ! ------------------------------------------------------------
      !    Search 4-Grid Points
      ! ------------------------------------------------------------      
      dist(1:num_grid_ave) = 1.D+10
      wk_maxp = num_grid_ave    
      DO jy=1,ny
        DO ix=1,nx
          CALL com_distll_1(flon(ix,jy),flat(ix,jy),&
                          & olon(inum) ,olat(inum) ,wk_dist)
          IF(wk_dist > max_dist) CYCLE
          IF(wk_dist < dist(wk_maxp)) THEN
            dist(wk_maxp) = wk_dist
            dist_min_x(inum, wk_maxp) = ix
            dist_min_y(inum, wk_maxp) = jy
            DO ip = 1, num_grid_ave
              IF(dist(ip) == maxval(dist(1:num_grid_ave))) THEN
                wk_maxp = ip
                EXIT
              END IF
            END DO
          END IF
        END DO
      END DO
      IF(detailout) WRITE(6,'(A,4(A,I4,A,I4,A))')  '  Intp Grids : ', &
        & '(', INT(dist_min_x(inum, 1)), ',', INT(dist_min_y(inum, 1)), ') ', &
        & '(', INT(dist_min_x(inum, 2)), ',', INT(dist_min_y(inum, 2)), ') ', &
        & '(', INT(dist_min_x(inum, 3)), ',', INT(dist_min_y(inum, 3)), ') ', &
        & '(', INT(dist_min_x(inum, 4)), ',', INT(dist_min_y(inum, 4)), ') '
      ! ------------------------------------------------------------
      !    Interpolation
      ! ------------------------------------------------------------ 
      dist(1:num_grid_ave) =  dist(1:num_grid_ave) * 1.0D-3
      sum_dist = dist(1) * dist(1) * dist(2) * dist(2) * dist(3) * dist(3)  &
             & + dist(2) * dist(2) * dist(3) * dist(3) * dist(4) * dist(4)  &
             & + dist(3) * dist(3) * dist(4) * dist(4) * dist(1) * dist(1)  &
             & + dist(4) * dist(4) * dist(1) * dist(1) * dist(2) * dist(2)
      ratio(1) = (dist(2)*dist(2)*dist(3)*dist(3)*dist(4)*dist(4))/sum_dist
      ratio(2) = (dist(3)*dist(3)*dist(4)*dist(4)*dist(1)*dist(1))/sum_dist
      ratio(3) = (dist(4)*dist(4)*dist(1)*dist(1)*dist(2)*dist(2))/sum_dist
      ratio(4) = (dist(1)*dist(1)*dist(2)*dist(2)*dist(3)*dist(3))/sum_dist
      IF(detailout) WRITE(6,'(2X,A,5F15.5)') 'ratio      :',ratio(1:4),SUM(ratio(1:4))
      oi(inum) = SUM(ratio(1:num_grid_ave) * dist_min_x(inum, 1:num_grid_ave))
      oj(inum) = SUM(ratio(1:num_grid_ave) * dist_min_y(inum, 1:num_grid_ave))
      IF(detailout) WRITE(6,'(2X,A,2F15.5)') 'position   :',oi(inum),oj(inum)
    END DO Obs_Loop_2
  END IF

  RETURN
END SUBROUTINE com_pos2ij
!-----------------------------------------------------------------------
! UTC to TAI93
!-----------------------------------------------------------------------
SUBROUTINE com_utc2tai(iy,im,id,ih,imin,sec,tai93)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: iy,im,id,ih,imin
  REAL(r_size),INTENT(IN) :: sec
  REAL(r_size),INTENT(OUT) :: tai93
  REAL(r_size),PARAMETER :: mins = 60.0d0
  REAL(r_size),PARAMETER :: hour = 60.0d0*mins
  REAL(r_size),PARAMETER :: day = 24.0d0*hour
  REAL(r_size),PARAMETER :: year = 365.0d0*day
  INTEGER,PARAMETER :: mdays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER :: days,i

  tai93 = REAL(iy-1993,r_size)*year + FLOOR(REAL(iy-1993)/4.0,r_size)*day
  days = id -1
  DO i=1,12
    IF(im > i) days = days + mdays(i)
  END DO
  IF(MOD(iy,4) == 0 .AND. im > 2) days = days + 1 !leap year
  tai93 = tai93 + REAL(days,r_size)*day + REAL(ih,r_size)*hour &
              & + REAL(imin,r_size)*mins + sec
  IF(iy > 1993 .OR. (iy==1993 .AND. im > 6)) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 1994 .OR. (iy==1994 .AND. im > 6)) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 1995) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 1997 .OR. (iy==1997 .AND. im > 6)) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 1998) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 2005) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 2008) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 2012 .OR. (iy==2012 .AND. im > 6)) tai93 = tai93 + 1.0d0 !leap second

  RETURN
END SUBROUTINE com_utc2tai
!-----------------------------------------------------------------------
! TAI93 to UTC
!-----------------------------------------------------------------------
SUBROUTINE com_tai2utc(tai93,iy,im,id,ih,imin,sec)
  IMPLICIT NONE
  INTEGER,PARAMETER :: n=8 ! number of leap seconds after Jan. 1, 1993
  INTEGER,PARAMETER :: leapsec(n) = (/  15638399,  47174400,  94608001,&
                                  &    141868802, 189302403, 410227204,&
                                  &    504921605, 615254406/)
  REAL(r_size),INTENT(IN) :: tai93
  INTEGER,INTENT(OUT) :: iy,im,id,ih,imin
  REAL(r_size),INTENT(OUT) :: sec
  REAL(r_size),PARAMETER :: mins = 60.0d0
  REAL(r_size),PARAMETER :: hour = 60.0d0*mins
  REAL(r_size),PARAMETER :: day = 24.0d0*hour
  REAL(r_size),PARAMETER :: year = 365.0d0*day
  INTEGER,PARAMETER :: mdays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  REAL(r_size) :: wk,tai
  INTEGER :: days,i,leap

  tai = tai93
  sec = 0.0d0
  DO i=1,n
    IF(FLOOR(tai93) == leapsec(i)+1) sec = 60.0d0 + tai93-FLOOR(tai93,r_size)
    IF(FLOOR(tai93) > leapsec(i)) tai = tai -1.0d0
  END DO
  iy = 1993 + FLOOR(tai /year)
  wk = tai - REAL(iy-1993,r_size)*year - FLOOR(REAL(iy-1993)/4.0,r_size)*day
  IF(wk < 0.0d0) THEN
    iy = iy -1
    wk = tai - REAL(iy-1993,r_size)*year - FLOOR(REAL(iy-1993)/4.0,r_size)*day
  END IF
  days = FLOOR(wk/day)
  wk = wk - REAL(days,r_size)*day
  im = 1
  DO i=1,12
    leap = 0
    IF(im == 2 .AND. MOD(iy,4)==0) leap=1
    IF(im == i .AND. days >= mdays(i)+leap) THEN
      im = im + 1
      days = days - mdays(i)-leap
    END IF
  END DO
  id = days +1

  ih = FLOOR(wk/hour)
  wk = wk - REAL(ih,r_size)*hour
  imin = FLOOR(wk/mins)
  IF(sec < 60.0d0) sec = wk - REAL(imin,r_size)*mins

  RETURN
END SUBROUTINE com_tai2utc
!-----------------------------------------------------------------------
! Date and time regulation
!-----------------------------------------------------------------------
SUBROUTINE com_datetime_reg(iy,im,id,ih,imin,isec)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: iy,im,id,ih,imin,isec
  INTEGER :: mdays

  DO WHILE(im <= 0)
    im = im + 12
    iy = iy - 1
  END DO
  DO WHILE(im > 12)
    im = im - 12
    iy = iy + 1
  END DO
  DO WHILE(isec < 0)
    isec = isec + 60
    imin = imin - 1
  END DO
  DO WHILE(isec >= 60)
    isec = isec - 60
    imin = imin + 1
  END DO
  DO WHILE(imin < 0)
    imin = imin + 60
    ih = ih - 1
  END DO
  DO WHILE(imin >= 60)
    imin = imin - 60
    ih = ih + 1
  END DO
  DO WHILE(ih < 0)
    ih = ih + 24
    id = id - 1
  END DO
  DO WHILE(ih >= 24)
    ih = ih - 24
    id = id + 1
  END DO
  DO WHILE(id <= 0)
    im = im - 1
    IF(im <= 0) THEN
      im = im + 12
      iy = iy - 1
    END IF
    CALL com_mdays(iy,im,mdays)
    id = id + mdays
  END DO
  CALL com_mdays(iy,im,mdays)
  DO WHILE(id > mdays)
    id = id - mdays
    im = im + 1
    IF(im > 12) THEN
      im = im - 12
      iy = iy + 1
    END IF
    CALL com_mdays(iy,im,mdays)
  END DO

  RETURN
END SUBROUTINE com_datetime_reg
!-----------------------------------------------------------------------
! Number of days of the month
!-----------------------------------------------------------------------
SUBROUTINE com_mdays(iy,im,mdays)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iy,im
  INTEGER, INTENT(OUT) :: mdays

  SELECT CASE(im)
  CASE(1,3,5,7,8,10,12)
    mdays = 31
  CASE(4,6,9,11)
    mdays = 30
  CASE(2)
    mdays = 28
    IF(MOD(iy,100) == 0) THEN
      IF(MOD(iy,400) == 0) THEN
        mdays = 29
      END IF
    ELSE IF(MOD(iy,4) == 0) THEN
      mdays = 29
    END IF
  CASE DEFAULT
    WRITE(6,'(A)') '[Error] com_mdays: invalid month.'
    STOP
  END SELECT
END SUBROUTINE com_mdays


!===============================================================================
recursive subroutine quick_sort_asnd(var,init,first,last)
implicit none
  integer :: first, last, i, j, it
  integer :: init(*)
  real(r_size) :: var(*) , x,t

  x = var( (first+last) / 2 )
  i = first  ;  j = last
  do
    do while (var(i) < x)  ;  i=i+1  ;  end do
    do while (x < var(j))  ;  j=j-1  ;  end do
    if (i >= j) exit
      t       = var(i)  ;  var(i)  = var(j)  ;  var(j)  = t
      it      = init(i) ;  init(i) = init(j) ;  init(j) = it
      i       = i + 1   ;  j       = j - 1
  end do
  if (first < i - 1 ) call quick_sort_asnd(var,init, first, i - 1)
  if (j + 1 < last  ) call quick_sort_asnd(var,init, j + 1, last )

  return
end subroutine quick_sort_asnd

!=============================================================================== SK from KK 20200408
!-----------------------------------------------------------------------
! Kullback–Leibler divergence (KLD)     ! KK 20160414
!-----------------------------------------------------------------------
subroutine com_kld(ndim, var1, kld)
  implicit none
  INTEGER, INTENT(IN)      :: ndim
  REAL(r_size),INTENT(IN)  :: var1(ndim)
  REAL(r_size),INTENT(OUT) :: kld
  REAL(r_size) :: hstg(ndim)
  REAL(r_size) :: work1(ndim)
  REAL(r_size) :: mean, var, std
  REAL(r_size) :: alpha, w
  REAL(r_size) :: threshold, x
  REAL(r_size) :: p, q, kld_tmp
  INTEGER :: ic, iclass, class_lmt
  INTEGER :: i

  kld = 0.d0
  work1(:) = var1(:)

  call qsort(ndim, work1)

  !--- calc. variance ---!
  mean = sum(work1(:))/dble(ndim)
  work1(:) = work1(:)-mean
  var = sum(work1(:)**2)/dble(ndim-1)
  if( var==0.0d0 )then
    kld= undef
    return
  endif
  std = sqrt(var)
  work1(:) = work1(:)/std

! do i = 1, ndim
!   print *, i, work1(i)
! end do

  !--- set histogram parameters ---!
  var = 1.d0
  std = 1.d0
  alpha = ndim**(1.d0/3.d0)
  w = 3.5d0*std/alpha

  !--- make histogram ---!
  iclass = 1
  do while(w*iclass+0.5d0*w < 10.d0*std)  !!!!  change 10.d0 - > 20.d0
!   print *, iclass, w*iclass+0.5d0*w, 10.d0*std
    iclass = iclass + 1
  end do
  class_lmt = iclass

  ic = 0
  iclass = 1
  hstg(:) = 0.d0
  i = 1
  do while(i <= ndim)
    threshold = -w*(class_lmt+0.5d0) + w*iclass
    if(work1(i) < threshold) then
      ic = ic + 1
!     print '(I5,I4,I6,2F11.5,F6.0)', ic, iclass, i, work1(i), threshold
    else
      hstg(iclass) = dble(ic)
!     print '(I5,I4,I6,2F11.5,F6.0)', ic, iclass, i, work1(i), threshold, hstg(iclass)
      iclass = iclass + 1
      ic = 0
      i = i - 1
    end if
    i = i + 1
  end do
  hstg(iclass) = dble(ic)
! print '(F8.0,I5)', sum(hstg(:)), iclass

  !--- normalize distribution ---!
  hstg(:) = hstg(:)/ndim

  !--- KL divergence ---!
  do iclass = 1, class_lmt*2+1
    x = w*iclass-w*(class_lmt+1)
    p = hstg(iclass)
    q = dexp(-x**2/(2*var))/dsqrt(2*pi*var)*w
    if(p == 0.d0) then
      kld_tmp = -p * log(q)
    else
      kld_tmp = p * (log(p) - log(q))
    end if
!   if(q == 0.d0) then
!     kld_tmp = p * log(p)
!   else
!     kld_tmp = p * (log(p) - log(q))
!   end if
    kld = kld + kld_tmp
!   print *, iclass, p, q
!   print '(I5,F7.0, F10.5)', iclass, hstg(iclass)*ndim, kld_tmp
  end do
  if(kld <= 0.d0) kld = 0.d0

  return
end subroutine com_KLD
!-----------------------------------------------------------------------
! KLD thresholds
!-----------------------------------------------------------------------
SUBROUTINE kld_threshold(nbv, kld_floor, kld_cap)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nbv
  REAL(r_size), INTENT(OUT) :: kld_floor, kld_cap
  REAL(r_size) :: kld_tmp

  kld_cap = 0.60d0
! kld_cap = 0.80d0

! select case(nbv)
!   case(:19)
!     kld_tmp = 0.25d0 - (0.25d0 - 0.17d0) / 20 * nbv
!   case(20:39)
!     kld_tmp = 0.17d0 - (0.17d0 - 0.13d0) / (40-20) * (nbv-20)
!   case(40:79)
!     kld_tmp = 0.13d0 - (0.13d0 - 0.09d0) / (80-40) * (nbv-40)
!   case(80:159)
!     kld_tmp = 0.09d0 - (0.09d0 - 0.06d0) / (160-80) * (nbv-80)
!   case(160:319)
!     kld_tmp = 0.06d0 - (0.06d0 - 0.05d0) / (320-160) * (nbv-160)
!   case(320:639)
!     kld_tmp = 0.05d0 - (0.05d0 - 0.03d0) / (640-320) * (nbv-320)
!   case(640:2559)
!     kld_tmp = 0.03d0 - (0.03d0 - 0.02d0) / (2560-640) * (nbv-640)
!   case(2560:)
!     kld_tmp = 0.02d0
! end select

! select case(nbv)
!   case(:19)
!     kld_tmp = 0.25d0 - (0.25d0 - 0.09d0) / 20 * nbv
!   case(20:39)
!     kld_tmp = 0.09d0 - (0.09d0 - 0.07d0) / (40-20) * (nbv-20)
!   case(40:79)
!     kld_tmp = 0.07d0 - (0.07d0 - 0.06d0) / (80-40) * (nbv-40)
!   case(80:159)
!     kld_tmp = 0.06d0 - (0.06d0 - 0.04d0) / (160-80) * (nbv-80)
!   case(160:319)
!     kld_tmp = 0.04d0 - (0.04d0 - 0.03d0) / (320-160) * (nbv-160)
!   case(320:1279)
!     kld_tmp = 0.03d0 - (0.03d0 - 0.02d0) / (1280-320) * (nbv-320)
!   case(1280:)
!     kld_tmp = 0.02d0
! end select

  !--- 75 percentile ---!
! select case(nbv)
!   case(:19)
!     kld_tmp = 0.250d0 - (0.25d0 - 0.086d0) / 20 * nbv
!   case(20:39)
!     kld_tmp = 0.086d0 - (0.086d0 - 0.066d0) / (40-20) * (nbv-20)
!   case(40:79)
!     kld_tmp = 0.066d0 - (0.066d0 - 0.049d0) / (80-40) * (nbv-40)
!   case(80:159)
!     kld_tmp = 0.049d0 - (0.049d0 - 0.034d0) / (160-80) * (nbv-80)
!   case(160:319)
!     kld_tmp = 0.034d0 - (0.034d0 - 0.023d0) / (320-160) * (nbv-160)
!   case(320:639)
!     kld_tmp = 0.023d0 - (0.023d0 - 0.015d0) / (640-320) * (nbv-320)
!   case(640:1279)
!     kld_tmp = 0.015d0 - (0.015d0 - 0.010d0) / (1280-640) * (nbv-640)
!   case(1280:)
!     kld_tmp = 0.010d0
! end select

  !--- 80 percentile ---!
! select case(nbv)
!   case(:19)
!     kld_tmp = 0.250d0 - (0.25d0 - 0.098d0) / 20 * nbv
!   case(20:39)
!     kld_tmp = 0.098d0 - (0.098d0 - 0.072d0) / (40-20) * (nbv-20)
!   case(40:79)
!     kld_tmp = 0.072d0 - (0.072d0 - 0.053d0) / (80-40) * (nbv-40)
!   case(80:159)
!     kld_tmp = 0.053d0 - (0.053d0 - 0.037d0) / (160-80) * (nbv-80)
!   case(160:319)
!     kld_tmp = 0.037d0 - (0.037d0 - 0.024d0) / (320-160) * (nbv-160)
!   case(320:639)
!     kld_tmp = 0.024d0 - (0.024d0 - 0.016d0) / (640-320) * (nbv-320)
!   case(640:1279)
!     kld_tmp = 0.016d0 - (0.016d0 - 0.011d0) / (1280-640) * (nbv-640)
!   case(1280:2559)
!     kld_tmp = 0.011d0 - (0.011d0 - 0.010d0) / (2560-1280) * (nbv-1280)
!   case(2560:)
!     kld_tmp = 0.010d0
! end select

  !--- 85 percentile ---!
! select case(nbv)
!   case(:19)
!     kld_tmp = 0.250d0 - (0.25d0 - 0.109d0) / 20 * nbv
!   case(20:39)
!     kld_tmp = 0.109d0 - (0.109d0 - 0.084d0) / (40-20) * (nbv-20)
!   case(40:79)
!     kld_tmp = 0.084d0 - (0.084d0 - 0.059d0) / (80-40) * (nbv-40)
!   case(80:159)
!     kld_tmp = 0.059d0 - (0.059d0 - 0.040d0) / (160-80) * (nbv-80)
!   case(160:319)
!     kld_tmp = 0.040d0 - (0.040d0 - 0.026d0) / (320-160) * (nbv-160)
!   case(320:639)
!     kld_tmp = 0.026d0 - (0.026d0 - 0.017d0) / (640-320) * (nbv-320)
!   case(640:1279)
!     kld_tmp = 0.017d0 - (0.017d0 - 0.011d0) / (1280-640) * (nbv-640)
!   case(1280:2559)
!     kld_tmp = 0.011d0 - (0.011d0 - 0.010d0) / (2560-1280) * (nbv-1280)
!   case(2560:)
!     kld_tmp = 0.010d0
! end select

  !--- 90 percentile ---!
  select case(nbv)
    case(:19)
      kld_tmp = 0.250d0 - (0.25d0 - 0.125d0) / 20 * nbv
    case(20:39)
      kld_tmp = 0.125d0 - (0.125d0 - 0.096d0) / (40-20) * (nbv-20)
    case(40:79)
      kld_tmp = 0.096d0 - (0.096d0 - 0.066d0) / (80-40) * (nbv-40)
    case(80:159)
      kld_tmp = 0.066d0 - (0.066d0 - 0.044d0) / (160-80) * (nbv-80)
    case(160:319)
      kld_tmp = 0.044d0 - (0.044d0 - 0.029d0) / (320-160) * (nbv-160)
    case(320:639)
      kld_tmp = 0.029d0 - (0.029d0 - 0.018d0) / (640-320) * (nbv-320)
    case(640:1279)
      kld_tmp = 0.018d0 - (0.018d0 - 0.012d0) / (1280-640) * (nbv-640)
    case(1280:2559)
      kld_tmp = 0.012d0 - (0.012d0 - 0.010d0) / (2560-1280) * (nbv-1280)
    case(2560:)
      kld_tmp = 0.010d0
  end select

  !--- 95 percentile ---!
! select case(nbv)
!   case(:19)
!     kld_tmp = 0.250d0 - (0.25d0 - 0.158d0) / 20 * nbv
!   case(20:39)
!     kld_tmp = 0.158d0 - (0.158d0 - 0.118d0) / (40-20) * (nbv-20)
!   case(40:79)
!     kld_tmp = 0.118d0 - (0.118d0 - 0.077d0) / (80-40) * (nbv-40)
!   case(80:159)
!     kld_tmp = 0.077d0 - (0.077d0 - 0.050d0) / (160-80) * (nbv-80)
!   case(160:319)
!     kld_tmp = 0.050d0 - (0.050d0 - 0.032d0) / (320-160) * (nbv-160)
!   case(320:639)
!     kld_tmp = 0.032d0 - (0.032d0 - 0.021d0) / (640-320) * (nbv-320)
!   case(640:1279)
!     kld_tmp = 0.021d0 - (0.021d0 - 0.013d0) / (1280-640) * (nbv-640)
!   case(1280:2559)
!     kld_tmp = 0.013d0 - (0.013d0 - 0.010d0) / (2560-1280) * (nbv-1280)
!   case(2560:)
!     kld_tmp = 0.010d0
! end select

! kld_floor = max(kld_tmp, 0.02d0)
  kld_floor = max(kld_tmp, 0.010d0)

  RETURN
END SUBROUTINE kld_threshold
!-----------------------------------------------------------------------
! Quick sort     ! KK 20171110
!-----------------------------------------------------------------------
recursive subroutine qsort(ndim, var)
  implicit none

  INTEGER, INTENT(IN)        :: ndim
  REAL(r_size),INTENT(INOUT) :: var(ndim)
  INTEGER :: i, imin, imax, is
  REAL(r_size) :: p, it

  if(ndim == 1) return

  p = var(1)
  imin = 1
  imax = ndim
  do while(imin < imax)
    do while(var(imin) < p .and. imin < imax)
      imin = imin + 1
    enddo
    do while(p <= var(imax) .and. imin < imax)
      imax = imax - 1
    enddo
    if(imin < imax)then
      it = var(imax)
      var(imax) = var(imin)
      var(imin) = it
    endif
  enddo
  is = imax - 1
  if(is == 0) is = 1
  call qsort(is,      var(1:is))
  call qsort(ndim-is, var(is+1:ndim))

  return
end subroutine qsort
!-----------------------------------------------------------------------
recursive subroutine iqsort(ndim, var)
  implicit none

  INTEGER, INTENT(IN)    :: ndim
  INTEGER, INTENT(INOUT) :: var(ndim)
  INTEGER :: i, imin, imax, is
  INTEGER :: p, it

  if(ndim == 1) return

  p = var(1)
  imin = 1
  imax = ndim
  do while(imin < imax)
    do while(var(imin) < p .and. imin < imax)
      imin = imin + 1
    enddo
    do while(p <= var(imax) .and. imin < imax)
      imax = imax - 1
    enddo
    if(imin < imax)then
      it = var(imax)
      var(imax) = var(imin)
      var(imin) = it
    endif
  enddo
  is = imax - 1
  if(is == 0) is = 1
  call iqsort(is,      var(1:is))
  call iqsort(ndim-is, var(is+1:ndim))

  return
end subroutine iqsort
!-----------------------------------------------------------------------
recursive subroutine qsort2(ndim,var,r)
  implicit none

  integer, intent(in) :: ndim
  real(8), intent(inout) :: var(ndim)
  integer, intent(inout) :: r(ndim)
  integer :: i,imin,imax,is
  real(8) :: p, it
  integer :: ir

  if(ndim == 1) return
  p = var(1)
  imin = 1
  imax = ndim
  do while(imin < imax)
    do while(var(imin) < p .and. imin < imax)
      imin = imin + 1
    enddo
    do while(p <= var(imax) .and. imin < imax)
      imax = imax - 1
    enddo
    if(imin < imax)then
      it = var(imax)
      ir = r(imax)
      var(imax) = var(imin)
      r(imax) = r(imin)
      var(imin) = it
      r(imin) = ir
    endif
  enddo
  is = imax - 1
  if(is == 0) is = 1
  call qsort2(is, var(1:is), r(1:is))
  call qsort2(ndim-is, var(is + 1:ndim), r(is + 1:ndim))

  return
end subroutine qsort2
!-----------------------------------------------------------------------
recursive subroutine iqsort2(ndim,var,r)
  implicit none

  integer, intent(in) :: ndim
  integer, intent(inout) :: var(ndim)
  integer, intent(inout) :: r(ndim)
  integer :: i,imin,imax,is
  integer :: p, it
  integer :: ir

  if(ndim == 1) return
  p = var(1)
  imin = 1
  imax = ndim
  do while(imin < imax)
    do while(var(imin) < p .and. imin < imax)
      imin = imin + 1
    enddo
    do while(p <= var(imax) .and. imin < imax)
      imax = imax - 1
    enddo
    if(imin < imax)then
      it = var(imax)
      ir = r(imax)
      var(imax) = var(imin)
      r(imax) = r(imin)
      var(imin) = it
      r(imin) = ir
    endif
  enddo
  is = imax - 1
  if(is == 0) is = 1
  call iqsort2(is, var(1:is), r(1:is))
  call iqsort2(ndim-is, var(is + 1:ndim), r(is + 1:ndim))

  return
end subroutine iqsort2
!-----------------------------------------------------------------------
recursive subroutine qsort3(ndim,var,r)
  implicit none

  integer, intent(in) :: ndim
  real(8), intent(inout) :: var(ndim)
  real(8), intent(inout) :: r(ndim)
  integer :: i,imin,imax,is
  real(8) :: p, it, ir

  if(ndim == 1) return
  p = var(1)
  imin = 1
  imax = ndim
  do while(imin < imax)
    do while(var(imin) < p .and. imin < imax)
      imin = imin + 1
    enddo
    do while(p <= var(imax) .and. imin < imax)
      imax = imax - 1
    enddo
    if(imin < imax)then
      it = var(imax)
      ir = r(imax)
      var(imax) = var(imin)
      r(imax) = r(imin)
      var(imin) = it
      r(imin) = ir
    endif
  enddo
  is = imax - 1
  if(is == 0) is = 1
  call qsort3(is, var(1:is), r(1:is))
  call qsort3(ndim-is, var(is + 1:ndim), r(is + 1:ndim))

  return
end subroutine qsort3
END MODULE common