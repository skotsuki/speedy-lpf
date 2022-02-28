program plot_station
  implicit none

  integer, parameter :: nlon = 96, nlat = 48
  real(4) :: temp = 10
  real(4) :: lon(nlon)
  real(4) :: lat(nlat)=(/-87.159,-83.479,-79.777,-76.070,-72.362,-68.652,-64.942,-61.232,-57.521,-53.810,-50.099,-46.389,-42.678,-38.967,-35.256,-31.545,-27.833,-24.122,-20.411,-16.700,-12.989, -9.278, -5.567, -1.856,  1.856,  5.567,  9.278,12.989, 16.700, 20.411, 24.122, 27.833, 31.545, 35.256, 38.967, 42.678, 46.389,50.099, 53.810, 57.521, 61.232, 64.942, 68.652, 72.362, 76.070, 79.777, 83.479, 87.159/)
  real(4) :: pi
  character(128) :: stid
  integer :: nflag, nlev
  real(4) :: tim
  integer :: i, j, nstation, cdummy
  integer :: ios
  integer :: ix, iy

  pi = 4.0*atan(1.0)
  do i = 1, nlon
    lon(i) = 0.0 + 3.75*(i-1)
  end do

  open(10, file='station.tbl.tmp')
  open(11, file='datafile.dat', form='unformatted', access='stream')
  open(12, file='obs_map.txt')

  nstation = 0
  read(10, '(A)') cdummy
  read(10, '(A)') cdummy
  do
    read(10, '(2I3)', iostat=ios) i, j
    if(ios /= 0) exit
    nstation = nstation + 1

    tim = 0.0
    nlev = 1
    nflag = 1
    write(stid, *) nstation
    write(11) trim(stid), j, i, tim, nlev, nflag
    write(11) temp
!   print *, i, j, nstation, trim(stid)
    write(12, *) lon(i), lat(j), "1"
    print *, nstation, i, j, lon(i), lat(j)
  end do

  nlev = 0
  write(11) trim(stid), j, i, tim, nlev, nflag

  close(10)
  close(11)
  close(12)


  stop
end program plot_station
