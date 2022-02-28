program prg_ensemble_mean
  implicit none

  integer, parameter :: nlev = 7, nlon = 96, nlat = 48
  integer, parameter :: nij0 = nlon*nlat
  integer, parameter :: nv3d = 4, nv2d = 2, pnv3d=5
  real(8) :: mean_3d(nlon, nlat, nlev, nv3d), mean_2d(nlon, nlat, nv2d)
  real(8) :: mean_p3d(nlon, nlat, nlev, pnv3d), mean_p2d(nlon, nlat, nv2d)
  real(8) :: sprd_3d(nlon, nlat, nlev, nv3d), sprd_2d(nlon, nlat, nv2d)
  real(8) :: sprd_p3d(nlon, nlat, nlev, pnv3d), sprd_p2d(nlon, nlat, nv2d)
  real(8), allocatable :: ens_3d(:,:,:,:,:), ens_2d(:,:,:,:)
  real(8), allocatable :: ens_p3d(:,:,:,:,:), ens_p2d(:,:,:,:)
  real(8) :: sum_3d(nlon, nlat, nlev, nv3d), sum_2d(nlon, nlat, nv2d)
  real(8) :: sum_p3d(nlon, nlat, nlev, pnv3d), sum_p2d(nlon, nlat, nv2d)
  integer :: member, l_member
  integer :: im
  character(10) :: ymdh
  character(256) :: ifile, ofile
  character(256) :: data_dir

  integer :: ifile_id
  integer :: ierr, nprocs, myrank

  include 'mpif.h'
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

  if(myrank == 0) print *, 'Compute Ensemble Mean & Spread with MPI'
  if(myrank == 0) print *, 'nprocs =', nprocs
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
! print *, 'myrank =', myrank

  ifile_id = 10 + myrank
 !open(ifile_id, file='info.txt')
 !read(ifile_id, '(A10)') ymdh
 !read(ifile_id, *) member
 !read(ifile_id, '(A256)') data_dir
 !close(ifile_id)
  if(myrank == 0) then
    open(10, file='info.txt')
    read(10, '(A10)') ymdh
    read(10, *) member
    read(10, '(A256)') data_dir
    close(10)
  end if
  CALL MPI_BCAST(ymdh, 10, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(member, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(data_dir, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
! if(myrank==1) print *, member, ymdh, data_dir


  if( mod(member, nprocs) == 0 ) then
    l_member = member/nprocs
  else
    l_member = idint(dint(dble(member)/dble(nprocs))) + 1
  end if

  allocate( ens_3d(nlon, nlat, nlev, nv3d, l_member) )
  allocate( ens_2d(nlon, nlat, nv2d, l_member) )
  allocate( ens_p3d(nlon, nlat, nlev, pnv3d, l_member) )
  allocate( ens_p2d(nlon, nlat, nv2d, l_member) )
  ens_3d(:,:,:,:,:) = 0.d0
  ens_2d(:,:,:,:) = 0.d0
  ens_p3d(:,:,:,:,:) = 0.d0
  ens_p2d(:,:,:,:) = 0.d0

  call read_data

  call calc_mean

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  call calc_sprd

  call output_data

  call MPI_Finalize(ierr)

contains
  SUBROUTINE read_grd4(filename,v3d,v2d)
    IMPLICIT NONE
    integer, parameter :: r_sngl = kind(0.0e0)
    CHARACTER(*),INTENT(IN) :: filename
    REAL(4),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
    REAL(4),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
    INTEGER :: iunit,iolen
    INTEGER :: i,j,k,n,irec

    iunit=11
    INQUIRE(IOLENGTH=iolen) iolen
    OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)
!   print *, r_sngl, iolend0

    irec=1
    DO n=1,nv3d
      DO k=1,nlev
        READ(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
        irec = irec + 1
      END DO
    END DO

    DO n=1,nv2d
      READ(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO

    CLOSE(iunit)

    RETURN
  END SUBROUTINE read_grd4
  !-----------------------------------------------------------------------!

  !-----------------------------------------------------------------------!
  SUBROUTINE write_grd4(filename,v3d,v2d)
    IMPLICIT NONE
    integer, parameter :: r_sngl = kind(0.0e0)
    CHARACTER(*),INTENT(IN) :: filename
    REAL(4),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
    REAL(4),INTENT(IN) :: v2d(nlon,nlat,nv2d)
    INTEGER :: iunit,iolen
    INTEGER :: i,j,k,n,irec

    iunit=11
    INQUIRE(IOLENGTH=iolen) iolen
    OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)
!   print *, r_sngl, iolen

    irec=1
    DO n=1,nv3d
      DO k=1,nlev
        write(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
        irec = irec + 1
      END DO
    END DO

    DO n=1,nv2d
      write(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO

    CLOSE(iunit)

    RETURN
  END SUBROUTINE write_grd4
  !-----------------------------------------------------------------------!

  !-----------------------------------------------------------------------!
  SUBROUTINE read_grd4_p(filename,v3d,v2d)
    IMPLICIT NONE
    integer, parameter :: r_sngl = kind(0.0e0)
    CHARACTER(*),INTENT(IN) :: filename
    REAL(4),INTENT(OUT) :: v3d(nlon,nlat,nlev,pnv3d)
    REAL(4),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
    INTEGER :: iunit,iolen
    INTEGER :: i,j,k,n,irec

    iunit=11
    INQUIRE(IOLENGTH=iolen) iolen
    OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)
!   print *, r_sngl, iolend0

    irec=1
    DO n=1,pnv3d
      DO k=1,nlev
        READ(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
        irec = irec + 1
      END DO
    END DO

    DO n=1,nv2d
      READ(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO

    CLOSE(iunit)

    RETURN
  END SUBROUTINE read_grd4_p
  !-----------------------------------------------------------------------!

  !-----------------------------------------------------------------------!
  SUBROUTINE write_grd4_p(filename,v3d,v2d)
    IMPLICIT NONE
    integer, parameter :: r_sngl = kind(0.0e0)
    CHARACTER(*),INTENT(IN) :: filename
    REAL(4),INTENT(IN) :: v3d(nlon,nlat,nlev,pnv3d)
    REAL(4),INTENT(IN) :: v2d(nlon,nlat,nv2d)
    INTEGER :: iunit,iolen
    INTEGER :: i,j,k,n,irec

    iunit=11
    INQUIRE(IOLENGTH=iolen) iolen
    OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=nij0*iolen)
!   print *, r_sngl, iolen

    irec=1
    DO n=1,pnv3d
      DO k=1,nlev
        write(iunit,REC=irec) ((v3d(i,j,k,n),i=1,nlon),j=1,nlat)
        irec = irec + 1
      END DO
    END DO

    DO n=1,nv2d
      write(iunit,REC=irec) ((v2d(i,j,n),i=1,nlon),j=1,nlat)
      irec = irec + 1
    END DO

    CLOSE(iunit)

    RETURN
  END SUBROUTINE write_grd4_p
  !-----------------------------------------------------------------------!

  !-----------------------------------------------------------------------!
  subroutine read_data
    implicit none
    character(128) :: ens_no
    integer :: im_tmp
    real(4) :: tmp_3d(nlon, nlat, nlev, nv3d), tmp_2d(nlon, nlat, nv2d)
    real(4) :: tmp_p3d(nlon, nlat, nlev, pnv3d), tmp_p2d(nlon, nlat, nv2d)

    !--- read model data ---!
    im_tmp = 1
    do im = 1 + myrank, member, nprocs
      write(ens_no, *) im
      if( im < 10 ) then
        ens_no = '00000'//trim(adjustl(ens_no))
      else if( im < 100 ) then
        ens_no = '0000'//trim(adjustl(ens_no))
      else if( im < 1000 ) then
        ens_no = '000'//trim(adjustl(ens_no))
      else if( im < 10000 ) then
        ens_no = '00'//trim(adjustl(ens_no))
      else if( im < 100000 ) then
        ens_no = '0'//trim(adjustl(ens_no))
      else
        ens_no = trim(adjustl(ens_no))
      end if
      ifile = trim(data_dir)//'/'//trim(ens_no)//'/'//ymdh//'.grd'
      call read_grd4(ifile, tmp_3d, tmp_2d)
      ens_3d(:,:,:,:,im_tmp) = dble(tmp_3d(:,:,:,:))
      ens_2d(:,:,:,im_tmp) = dble(tmp_2d(:,:,:))
      im_tmp = im_tmp + 1
    end do

    !--- read p data ---!
    im_tmp = 1
    do im = 1 + myrank, member, nprocs
      write(ens_no, *) im
      if( im < 10 ) then
        ens_no = '00000'//trim(adjustl(ens_no))
      else if( im < 100 ) then
        ens_no = '0000'//trim(adjustl(ens_no))
      else if( im < 1000 ) then
        ens_no = '000'//trim(adjustl(ens_no))
      else if( im < 10000 ) then
        ens_no = '00'//trim(adjustl(ens_no))
      else if( im < 100000 ) then
        ens_no = '0'//trim(adjustl(ens_no))
      else
        ens_no = trim(adjustl(ens_no))
      end if
      ifile = trim(data_dir)//'/'//trim(ens_no)//'/'//ymdh//'_p.grd'
      call read_grd4_p(ifile, tmp_p3d, tmp_p2d)
!     print *, myrank, im
      ens_p3d(:,:,:,:,im_tmp) = dble(tmp_p3d(:,:,:,:))
      ens_p2d(:,:,:,im_tmp) = dble(tmp_p2d(:,:,:))
      im_tmp = im_tmp + 1
    end do

!   im_tmp = 1
!   do im = 1 + myrank, member, nprocs
!     print *, myrank, im, im_tmp, ens_p3d(1,1,4,5,im_tmp)
!     im_tmp = im_tmp + 1
!   end do

    return
  end subroutine read_data
  !-----------------------------------------------------------------------!

  !-----------------------------------------------------------------------!
  subroutine calc_mean
    implicit none
    integer :: count, datatype, op, comm, ierr

    mean_3d(:,:,:,:) = 0.d0
    mean_2d(:,:,:) = 0.d0
    mean_p3d(:,:,:,:) = 0.d0
    mean_p2d(:,:,:) = 0.d0
    sum_3d(:,:,:,:) = 0.d0
    sum_2d(:,:,:) = 0.d0
    sum_p3d(:,:,:,:) = 0.d0
    sum_p2d(:,:,:) = 0.d0

    do im = 1, l_member
      sum_3d(:,:,:,:) = sum_3d(:,:,:,:) + ens_3d(:,:,:,:,im)
      sum_2d(:,:,:) = sum_2d(:,:,:) + ens_2d(:,:,:,im)
      sum_p3d(:,:,:,:) = sum_p3d(:,:,:,:) + ens_p3d(:,:,:,:,im)
      sum_p2d(:,:,:) = sum_p2d(:,:,:) + ens_p2d(:,:,:,im)
    end do

    count = nlon * nlat * nlev * nv3d
    call MPI_ALLREDUCE(sum_3d, mean_3d, count, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    MPI_COMM_WORLD, ierr)
    count = nlon * nlat * nv2d
    call MPI_ALLREDUCE(sum_2d, mean_2d, count, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    MPI_COMM_WORLD, ierr)
    count = nlon * nlat * nlev * pnv3d
    call MPI_ALLREDUCE(sum_p3d, mean_p3d, count, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    MPI_COMM_WORLD, ierr)
    count = nlon * nlat * nv2d
    call MPI_ALLREDUCE(sum_p2d, mean_p2d, count, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    MPI_COMM_WORLD, ierr)

    mean_3d(:,:,:,:) = mean_3d(:,:,:,:)/dble(member)
    mean_2d(:,:,:) = mean_2d(:,:,:)/dble(member)
    mean_p3d(:,:,:,:) = mean_p3d(:,:,:,:)/dble(member)
    mean_p2d(:,:,:) = mean_p2d(:,:,:)/dble(member)

!   print *, myrank, mean_p3d(1,1,4,5), mean_p3d(1,1,4,2)
!   do im = 1, l_member
!     print *, myrank, sum_p3d(1,1,4,5), ens_p3d(1,1,4,5,im)
!   end do

    return
  end subroutine calc_mean
  !-----------------------------------------------------------------------!

  !-----------------------------------------------------------------------!
  subroutine calc_sprd
    implicit none
    integer :: count, datatype, op, comm, ierr

    sprd_3d(:,:,:,:) = 0.d0
    sprd_2d(:,:,:) = 0.d0
    sprd_p3d(:,:,:,:) = 0.d0
    sprd_p2d(:,:,:) = 0.d0
    sum_3d(:,:,:,:) = 0.d0
    sum_2d(:,:,:) = 0.d0
    sum_p3d(:,:,:,:) = 0.d0
    sum_p2d(:,:,:) = 0.d0

    do im = 1, l_member
      if( sum(ens_3d(:,:,1,1,im)) /= 0.d0 ) then
        sum_3d(:,:,:,:) &
          = sum_3d(:,:,:,:) + (ens_3d(:,:,:,:,im) - mean_3d(:,:,:,:))**2
        sum_2d(:,:,:) &
          = sum_2d(:,:,:) + (ens_2d(:,:,:,im) - mean_2d(:,:,:))**2
        sum_p3d(:,:,:,:) &
          = sum_p3d(:,:,:,:) + (ens_p3d(:,:,:,:,im) - mean_p3d(:,:,:,:))**2
        sum_p2d(:,:,:) &
          = sum_p2d(:,:,:) + (ens_p2d(:,:,:,im) - mean_p2d(:,:,:))**2
      end if
    end do
    
    count = nlon * nlat * nlev * nv3d
    call MPI_REDUCE(sum_3d, sprd_3d, count, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                    MPI_COMM_WORLD, ierr)
    count = nlon * nlat * nv2d
    call MPI_REDUCE(sum_2d, sprd_2d, count, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                    MPI_COMM_WORLD, ierr)
    count = nlon * nlat * nlev * pnv3d
    call MPI_REDUCE(sum_p3d, sprd_p3d, count, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                    MPI_COMM_WORLD, ierr)
    count = nlon * nlat * nv2d
    call MPI_REDUCE(sum_p2d, sprd_p2d, count, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                    MPI_COMM_WORLD, ierr)

    if( myrank == 0 ) then
      sprd_3d(:,:,:,:) = dsqrt(sprd_3d(:,:,:,:)/dble(member-1))
      sprd_2d(:,:,:) = dsqrt(sprd_2d(:,:,:)/dble(member-1))
      sprd_p3d(:,:,:,:) = dsqrt(sprd_p3d(:,:,:,:)/dble(member-1))
      sprd_p2d(:,:,:) = dsqrt(sprd_p2d(:,:,:)/dble(member-1))
    end if

    return
  end subroutine calc_sprd
  !-----------------------------------------------------------------------!

  !-----------------------------------------------------------------------!
  subroutine output_data
    implicit none
    real(4) :: tmp_3d(nlon, nlat, nlev, nv3d), tmp_2d(nlon, nlat, nv2d)
    real(4) :: tmp_p3d(nlon, nlat, nlev, pnv3d), tmp_p2d(nlon, nlat, nv2d)

!   if(myrank==0) then
!     print *, sum(mean_3d(:,:,1,1))/(nlon*nlat), sum(mean_3d(:,:,1,4))/(nlon*nlat)
!     print *, sum(sprd_3d(:,:,1,1))/(nlon*nlat), sum(sprd_3d(:,:,1,4))/(nlon*nlat)
!   end if

    if( myrank == 0 ) then
      !--- write mixed analysis mean data (model) ---!
      ifile = trim(data_dir)//'/mean/'//trim(ymdh)//'.grd'
      tmp_3d(:,:,:,:) = real((mean_3d(:,:,:,:)))
      tmp_2d(:,:,:) = real(mean_2d(:,:,:))
      call write_grd4(ifile, tmp_3d, tmp_2d)

      !--- write mixed analysis spread data (model) ---!
      ifile = trim(data_dir)//'/sprd/'//trim(ymdh)//'.grd'
      tmp_3d(:,:,:,:) = real(sprd_3d(:,:,:,:))
      tmp_2d(:,:,:) = real(sprd_2d(:,:,:))
      call write_grd4(ifile, tmp_3d, tmp_2d)

      !--- write mixed analysis mean data (p) ---!
      ifile = trim(data_dir)//'/mean/'//trim(ymdh)//'_p.grd'
      tmp_p3d(:,:,:,:) = real(mean_p3d(:,:,:,:))
      tmp_p2d(:,:,:) = real(mean_p2d(:,:,:))
      call write_grd4_p(ifile, tmp_p3d, tmp_p2d)

      !--- write mixed analysis spread data (p) ---!
      ifile = trim(data_dir)//'/sprd/'//trim(ymdh)//'_p.grd'
      tmp_p3d(:,:,:,:) = real(sprd_p3d(:,:,:,:))
      tmp_p2d(:,:,:) = real(sprd_p2d(:,:,:))
      call write_grd4_p(ifile, tmp_p3d, tmp_p2d)
    end if

    return
  end subroutine output_data
  !-----------------------------------------------------------------------!
end program prg_ensemble_mean
