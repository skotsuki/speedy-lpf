program main
  implicit none

  !===> setting of speedy
  integer, parameter :: nlon  = 96
  integer, parameter :: nlat  = 48
  integer, parameter :: nlev  =  7
  integer, parameter :: nhour =  6
  real(8), parameter :: dx    = 3.75d0

  !===> setting of constant parameter
  real(4), parameter :: undef4 = -9.99d33
  real(8), parameter :: undef8 = -9.99d33
  real(8), parameter :: pi     =  3.1415926535d0
  real(8), parameter :: pi_180 = pi / 180.0d0
  REAL(8), PARAMETER :: lfactor     = 2.0d0*dsqrt(10.0d0/3.0d0)

  !===>speedy.cnf
  integer,        save :: adate, sdate, edate, pdate, qdate
  integer,        save :: nbv, nrunme
  character(255), save :: natdir, spndir, expdir, outtime, outtave, expname, obsname, cneff, crtps

  !===> save variables
  integer, save :: nsamp, nrec, nmax

  !====> global variables
  integer :: idate, m, mm, i, j, ii, jj, k, kk, irec, icon, igrd, timer
  real(8) :: asamp, psamp, www
  real(8) :: clon(nlon), clat(nlat)

  character(255) :: truname, emename, espname, pefname, outfile
  character(2)   :: CRGN
  logical :: ex

  integer, parameter :: tmax = 99999
  real(8) :: epeff, epeff_tave, epeff_t(tmax)
  real(8) :: ermse, ermse_tave, ermse_t(tmax)
  real(8) :: esprd, esprd_tave, esprd_t(tmax) 

  real(4) :: tru(nlon,nlat), eme(nlon,nlat), esp(nlon,nlat), erm(nlon,nlat)
  real(4) :: pef(nlon,nlat)

  real(4) :: erm_a, esp_a


  real(4) :: anal(nlon,nlat), true(nlon,nlat), rmse(nlon,nlat), sprd(nlon,nlat), diff(nlon,nlat)

  real(4) :: tmp4(nlon,nlat)
  real(8) :: tmp8(nlon,nlat)


  namelist / speedy_param / adate, sdate, edate, pdate, qdate, nbv, nrunme, &
                            natdir, spndir, expdir, outtime, outtave, expname, obsname, cneff, crtps
  !================================================================================================== (1) PREPARATION
  open(1,file="speedy.cnf")
  	read(1,nml=speedy_param)
  close(1)

  !=====> (1) preparation
  !print '(2a)', "  Calc... Standard  RMSE"
  do i=1,nlon  ;  clon(i) = 0.0d0 + dble(i-1) * dx   ;  end do
  call get_clat( nlat, clat )
  !==================================================================================================

  print '(3a)',">>>>>================================================   calc... RMSE&SPRD  ===============================================<<<<< "
  !!write(outfile,'(6a)') trim(outdir),"/time_",trim(obsname),"_",trim(expname),".txt"
  !!open(1,file=trim(outfile),form="formatted",action="write")
  open(1,file=trim(outtime),form="formatted",action="write")

  !=====> (2) get nsamp & RMSE
  idate = adate ! initial time
  nsamp = 0
  timer = 0
  erm_a = 0.0d0  ;  esp_a = 0.0d0

  rmse  = 0.0d0
  sprd  = 0.0d0
  ermse_t(:) = undef8  ;  esprd_t(:) = undef8
  epeff_t(:) = undef8
  do
    if( mod( idate,10000 )==100 ) &
      print '(5a,i10)',"      calc... RMSE&SPRD  ===>> ", trim(obsname)," ",trim(expdir)," :: ",idate

    timer = timer + 1
      write(truname,'(1a,i10.10,a)') trim(natdir),                 idate,".grd"
    if( idate < sdate )then
      write(emename,'(2a,i10.10,a)') trim(spndir),"/gues/mean/",   idate,".grd"
      write(espname,'(2a,i10.10,a)') trim(spndir),"/gues/sprd/",   idate,".grd"
    else
      write(emename,'(2a,i10.10,a)') trim(expdir),"/gues/mean/",   idate,".grd"
      write(espname,'(2a,i10.10,a)') trim(expdir),"/gues/sprd/",   idate,".grd"
      write(pefname,'(2a,i10.10,a)') trim(expdir),"/peff_lpf/",    idate,".grd"
    endif
    !!!print '(i10,x,a)', idate, trim(ename)

    ! U(7), V(7), T(7), P(7), Ps(2), Rain(2)
    k = 18 ! 4th level of T
    !!!k =  4 ! 4th level of U (to compare with Kondo and Miyoshi 2016)

    open(10,file=trim(truname),form="unformatted",access="direct",recl=nlon*nlat*4,action="read")
      read(10,rec=k) tru(:,:)
    close(10)

    inquire(FILE=trim(emename),EXIST=ex)

    if( ex .or. idate<=(sdate+700) ) then ! at least one week
      open(11,file=trim(emename),form="unformatted",access="direct",recl=nlon*nlat*4,action="read")
      open(12,file=trim(espname),form="unformatted",access="direct",recl=nlon*nlat*4,action="read")
        read(11,rec=k) eme(:,:)
        read(12,rec=k) esp(:,:)
      close(11) ; close(12)  
    else
      eme(:,:) = undef4
      esp(:,:) = undef4
      print '(2a)', " maybe because of divergence, cannot find file :: ", trim(emename)
    endif


    if( idate >= sdate ) then
      inquire(FILE=trim(pefname),EXIST=ex)
      if( ex )then
        open(20,file=trim(pefname),form="unformatted",access="direct",recl=nlon*nlat*4,action="read")
          read(20,rec=k) pef(:,:)
        close(20)
      else
        pef(:,:) = undef4
      end if
    else
      pef(:,:) = dble(nbv)
    end if

    !===> for all time (global-mean rmse)   
    ermse = 0.0d0 ; esprd = 0.0d0
    epeff = 0.0d0
    asamp = 0.0d0 ; psamp = 0.0d0
    do j=1,nlat
      www = dcos( clat(j)*pi_180 )
      do i=1,nlon
          asamp = asamp + www
          ermse = ermse + www * ( eme(i,j) - tru(i,j) )**2.0d0 
          esprd = esprd + www * ( esp(i,j)            )**2.0d0
        if( 0.0d0<=pef(i,j) .and. pef(i,j)<=dble(nbv) )then
          psamp = psamp + www
          epeff = epeff + www *   pef(i,j)
        end if
      end do
    end do

    ermse = dsqrt( ermse / asamp )  ;  esprd = dsqrt( esprd / asamp )    
    epeff =        epeff / psamp
    ermse_t(timer) = ermse          ;  esprd_t(timer) = esprd
    epeff_t(timer) = epeff

    ermse_tave = undef8 ; if( timer >=nrunme ) ermse_tave = sum( ermse_t(timer-nrunme+1:timer) ) / dble( nrunme )
    esprd_tave = undef8 ; if( timer >=nrunme ) esprd_tave = sum( esprd_t(timer-nrunme+1:timer) ) / dble( nrunme ) 
    epeff_tave = undef8 ; if( timer >=nrunme ) epeff_tave = sum( epeff_t(timer-nrunme+1:timer) ) / dble( nrunme )

    if( timer>= nrunme ) &
    write(1,'(2i12,f9.3,4f8.4,2f8.2)')                                     &
                                idate, timer, real(timer)/4.0d0,           & 
                                ermse, esprd, ermse_tave, esprd_tave, epeff, epeff_tave
    !write(6,'(2i12,f9.3,4f8.4,2f8.2)')                                     &
    !                            idate, timer, real(timer)/4.0d0,           & 
    !                            ermse, esprd, ermse_tave, esprd_tave, epeff, epeff_tave

    !===> for time-mean rmse
    if( pdate<=idate .and. idate<=qdate )then
      nsamp = nsamp + 1
      erm_a = erm_a + ermse
      esp_a = esp_a + esprd
    endif

    call update_date(idate,nhour)
    if ( nsamp .ge. 10000 ) then
      write(6,*) 'nsamp is too large'
      stop
    end if
    if( idate > edate  ) goto 10
  end do
10 continue
  close(1)

  erm_a = erm_a / dble( nsamp )
  esp_a = esp_a / dble( nsamp )
  if( .not. (0.0d0<abs(erm_a) .and. abs(erm_a)<999.9d0) ) erm_a = 999.9d0  
  if( .not. (0.0d0<abs(esp_a) .and. abs(esp_a)<999.9d0) ) esp_a = 999.9d0
!  write(outfile,'(6a,i10.10,a,i10.10,a)') trim(outdir),"/tave_",trim(obsname),"_",trim(expname),"_SMP",sdate,"-",edate,".txt"
!  open(1,file=trim(outfile),form="formatted",action="write")
  open(1,file=trim(outtave),form="formatted",action="write")
    write(1,'(a,x,a,x,i12,4f10.4)') trim(cneff), trim(crtps), nsamp, erm_a, esp_a
  close(1)


1111 continue

end
!=========================================================
subroutine calc_dist_deg(xlon1,xlon2,xlat1,xlat2,dist)
!=========================================================
  implicit none
  double precision xlon1,xlon2,xlat1,xlat2
  double precision lon1,lon2,lat1,lat2,cosd,dist

  REAL(8),PARAMETER :: pi=3.1415926535d0
  REAL(8),PARAMETER :: re=6371.3d3
  REAL(8),PARAMETER :: r180=1.0d0/180.0d0

  lon1 = xlon1 * pi * r180 ! [rad]
  lon2 = xlon2 * pi * r180 ! [rad]
  lat1 = xlat1 * pi * r180 ! [rad]
  lat2 = xlat2 * pi * r180 ! [rad]

  cosd = SIN(lat1)*SIN(lat2) + COS(lat1)*COS(lat2)*COS(lon2-lon1)
  cosd = MIN( 1.d0,cosd)
  cosd = MAX(-1.d0,cosd)
  dist = ACOS( cosd ) * re

end
!=========================================================
subroutine update_date(date,nhour)
!=========================================================
  implicit none
  integer date,nhour,date1
  integer year,mon,day,hour,mondays

  date1 = date
  year = int( date / 1000000 )  ;  date = date - year * 1000000
  mon  = int( date / 10000 )    ;  date = date - mon  * 10000
  day  = int( date / 100 )      ;  hour = date - day * 100

  hour = hour + nhour

10 continue
  if( hour.ge.24 ) then
    hour = hour - 24
    day  = day + 1
    call get_mondays(year,mon,mondays)
    if ( day.gt.mondays ) then
      day = 1
      mon = mon + 1
      if( mon.ge.13 ) then
        mon = 1
        year = year + 1
      end if ! mon
    end if ! day
  end if ! hour
  if( hour.ge.24 ) goto 10

  date = year*1000000 + mon*10000 + day*100 + hour
end subroutine update_date

!=========================================================
subroutine get_mondays(year,mon,mondays)
!=========================================================
  implicit none
  integer year,mon,mondays

  mondays=31
  if( mon ==  4 ) mondays = 30
  if( mon ==  6 ) mondays = 30
  if( mon ==  9 ) mondays = 30
  if( mon == 11 ) mondays = 30
  if( mon ==  2 ) then
    mondays = 28
    if( mod(year,4)==0 ) mondays = 29
  end if
end subroutine get_mondays
!=========================================================
subroutine get_clat(nlat,clat)
!=========================================================
implicit none
  integer, intent(in)  :: nlat
  real(8), intent(out) :: clat(nlat)

  clat(01)=-87.159  ;  clat(02)=-83.479  ;  clat(03)=-79.777  ;  clat(04)=-76.070  ;  clat(05)=-72.362  ;  clat(06)=-68.652  
  clat(07)=-64.942  ;  clat(08)=-61.232  ;  clat(09)=-57.521  ;  clat(10)=-53.810  ;  clat(11)=-50.099  ;  clat(12)=-46.389
  clat(13)=-42.678  ;  clat(14)=-38.967  ;  clat(15)=-35.256  ;  clat(16)=-31.545  ;  clat(17)=-27.833  ;  clat(18)=-24.122
  clat(19)=-20.411  ;  clat(20)=-16.700  ;  clat(21)=-12.989  ;  clat(22)=- 9.278  ;  clat(23)=- 5.567  ;  clat(24)=- 1.856
  clat(25)=  1.856  ;  clat(26)=  5.567  ;  clat(27)=  9.278  ;  clat(28)= 12.989  ;  clat(29)= 16.700  ;  clat(30)= 20.411
  clat(31)= 24.122  ;  clat(32)= 27.833  ;  clat(33)= 31.545  ;  clat(34)= 35.256  ;  clat(35)= 38.967  ;  clat(36)= 42.678
  clat(37)= 46.389  ;  clat(38)= 50.099  ;  clat(39)= 53.810  ;  clat(40)= 57.521  ;  clat(41)= 61.232  ;  clat(42)= 64.942
  clat(43)= 68.652  ;  clat(44)= 72.362  ;  clat(45)= 76.070  ;  clat(46)= 79.777  ;  clat(47)= 83.479  ;  clat(48)= 87.159
end subroutine get_clat
