program main
  implicit none

  !===> setting of speedy
  integer, parameter :: nbv    = 40
  integer, parameter :: ncase  = 1     ! the product (1000)
  integer, parameter :: ntime  = 100   ! the product (10000)

  integer :: init(nbv,nbv)
  real(8) :: rand(nbv), pfwgh(nbv), acc(nbv)
  real(8) :: pmat_answ(nbv,nbv), pmat_time(nbv,nbv,ntime), pmat_kkkk(nbv,nbv,ntime), tmp8(nbv,nbv)
  real(8) ::                                               pvec_kkkk(nbv,    ntime) 
  real(8) :: aerr(ncase,ntime) , perr(ncase,ntime)
  real(8) :: pmat_one(nbv,nbv), pmat(nbv,nbv), pmat_suu(nbv,nbv), pmat_mlt(nbv,nbv)
  
  !===> merginal particle filter
  integer, parameter :: nsu = 100

  !===> global variables
  integer :: i, j, k, l, itime, icase
  real(8) :: p_ave, p_std, m_ave, m_std, timer, peff
  character(255) :: fname

  integer t0, t1, t_rate, t_max, diff
  integer :: a_sample, b_sample, c_sample, d_sample
!  real(8), parameter :: a_factor =  0.088d0, b_factor = -0.088d0 
  real(8), parameter :: a_factor =  0.04d0, b_factor = -0.04d0
!  real(8), parameter :: a_factor =  0.040d0, b_factor = -0.040d0
!  real(8), parameter :: a_factor =  0.010d0, b_factor = -0.010d0

  real(8) :: c_factor, d_factor, sss1, sss2, ppp, qqq
!-----------------------------------------------------------------------
  ! MPF simple check
  a_sample = 100
  b_sample =  98

  sss1 = 0.0d0
  sss2 = 0.0d0
  do i=1,a_sample
    sss1 = sss1 + a_factor
    sss2 = sss2 + a_factor**2.0d0  ;  end do
  do i=1,b_sample
    sss1 = sss1 + b_factor
    sss2 = sss2 + b_factor**2.0d0  ;  end do

  ppp = 1.0d0 - sss1
  qqq = 1.0d0 - sss2
  c_factor = 0.5d0 * ( ppp - dsqrt(  ppp**2.0d0 - 2.0d0*(ppp**2.0d0 - qqq) ) )
  d_factor = 0.5d0 * ( ppp + dsqrt(  ppp**2.0d0 - 2.0d0*(ppp**2.0d0 - qqq) ) )


  sss1 = sss1 + c_factor        + d_factor 
  sss2 = sss2 + c_factor**2.0d0 + d_factor**2.0d0

  print '(4f)', a_factor, b_factor, c_factor, d_factor
  print '(4f)', sss1, sss2, ppp, qqq
  !!stop
!-----------------------------------------------------------------------
  aerr(:,:) = 0.0d0

  call system_clock(t0)
  DO ICASE=1,NCASE
    call random_number(rand)
    pfwgh(:) = 1.0d0          / dble( nbv )    ! (1) Peff = large
    pfwgh(:) = rand(:)        / sum( rand(:) ) ! (2) Peff = intermidiate
    rand(:)  = rand(:)**2.0d0                  ! (3) Peff = small
    pfwgh(:) = rand(:)        / sum( rand(:) ) ! (3) Peff = small

    peff     = 0.0d0
    do j=1,nbv
      tmp8(:,j) = pfwgh(:)
      peff      = peff + pfwgh(j)**2.0d0
    end do
    peff = 1.0d0 / peff

    !===> gen weight
    acc(:)   = 0.0d0
    acc(1)   = pfwgh(1)
    DO j=2,nbv  ;  acc(j) = acc(j-1) + pfwgh(j)  ;  END DO

    pmat(:,:) = 0.0d0
    DO j=1,a_sample ! MONTE-Carlo Resampling M times
      CALL get_resampling_mtx('MR','OF',nbv,acc,pmat_one)
      pmat(:,:) = pmat(:,:) + pmat_one(:,:)*a_factor ; END DO
    DO j=1,b_sample ! MONTE-Carlo Resampling M times
      CALL get_resampling_mtx('MR','OF',nbv,acc,pmat_one)
      pmat(:,:) = pmat(:,:) + pmat_one(:,:)*b_factor ; END DO

      CALL get_resampling_mtx('MR','OF',nbv,acc,pmat_one)
      pmat(:,:) = pmat(:,:) + pmat_one(:,:)*c_factor
      pmat_mlt(:,:) = pmat_one(:,:)

      !!!CALL get_resampling_mtx('SU','OF',nbv,acc,pmat_one)
      !!!!CALL get_resampling_mtx('MR','OF',nbv,acc,pmat_one)
    pmat_suu(:,:) = 0.0d0
    DO j=1,nsu
      CALL get_resampling_mtx('SU','OF',nbv,acc,pmat_one)
      pmat_suu(:,:) = pmat_suu(:,:) + pmat_one(:,:) / dble(nsu) ; END DO

      pmat(:,:) = pmat(:,:) + pmat_suu(:,:)*d_factor    


    !IF( icase==1 )THEN
      write(fname,'(a,i6.6,a,i8.8,a,i8.8,a)')"./out_900/sampled-error_M",nbv,"NTIME",ntime,"NCASE",ncase,".bin"
      open(1,file=trim(fname),form="unformatted",access="direct",recl=nbv*nbv*4,action="write")
        write(1,rec=1) (( sngl( tmp8(i,j)     ),j=1,nbv),i=1,nbv)
        write(1,rec=2) (( sngl( pmat_mlt(i,j) ),j=1,nbv),i=1,nbv)
        write(1,rec=3) (( sngl( pmat_one(i,j) ),j=1,nbv),i=1,nbv)
        write(1,rec=4) (( sngl( pmat_suu(i,j) ),j=1,nbv),i=1,nbv)
        write(1,rec=5) (( sngl( pmat(i,j)     ),j=1,nbv),i=1,nbv)
      close(1)
    !ENDIF
  END DO

    print '(a)', "========================ave time ============================="
    do i=1,nbv
      print '(I3,a,3f8.3,a,3f8.3)', i, &
        "   || one-su, mlt-su, ans :: ",sum(pmat_one(i,:)), sum(pmat_suu(i,:)), sum(pmat(i,:)), &
        "   || wgh, acc, peff :: ",pfwgh(i)*dble(nbv), acc(i), peff
    end do

!NONEED!  !==> statistics
!NONEED!  write(fname,'(a,i6.6,a,i8.8,a,i8.8,a)')"./out_900/sampled-error_M",nbv,"NTIME",ntime,"NCASE",ncase,".txt"
!NONEED!  open(1,file=trim(fname),form="formatted",action="write")
!NONEED!  DO itime=1,ntime
!NONEED!    p_ave = sum( perr(1:ncase,itime) ) / dble( NCASE )
!NONEED!    p_std = dsqrt( sum( (perr(1:ncase,itime)-p_ave)**2.0d0 ) / dble( NCASE ) )
!NONEED!    m_ave = sum( aerr(1:ncase,itime) ) / dble( NCASE )
!NONEED!    m_std = dsqrt( sum( (aerr(1:ncase,itime)-m_ave)**2.0d0 ) / dble( NCASE ) )
!NONEED!
!NONEED!    write(1,'(2i,6f15.8,4x,6f15.8)') ncase, itime,                                                          & 
!NONEED!      p_ave, p_std, p_ave-p_std, p_ave+p_std, minval( perr(1:ncase,itime) ), maxval( perr(1:ncase,itime) ), &
!NONEED!      m_ave, m_std, m_ave-m_std, m_ave+m_std, minval( aerr(1:ncase,itime) ), maxval( aerr(1:ncase,itime) )
!NONEED!    !!write(6,'(2i,4f15.8)') ncase, itime, ave, std, ave-std, ave+std
!NONEED!  END DO
!NONEED!  close(1)

  print '(a)', "program end from prg_demo-mptrns.f90"
  stop
end

!-----------------------------------------------------------------------
SUBROUTINE get_resampling_mtx(CC,DG,nbv,acc,pmat)
!-----------------------------------------------------------------------
IMPLICIT NONE
  INTEGER, INTENT(in)  :: nbv
  CHARACTER(2), INTENT(in)  :: CC, DG
  REAL(8), INTENT(in)  :: acc(nbv)
  REAL(8), INTENT(out) :: pmat(nbv,nbv)

  INTEGER :: i, j, k, init(nbv), inum(nbv)
  REAL(8) :: rand(nbv), temp

  DO j=1,nbv ; init(j) = j  ; END DO
  IF( CC == 'SU' )THEN
    !CALL com_rand_seed(nbv,0,rand) ! [0-1]
    CALL random_number(rand)
    temp = rand(1)
    DO i=1,nbv
      rand(i) = ( dble(i) +  temp - 1.0d0 ) / dble(nbv)
    END DO
  ELSE IF( CC == 'MR' )THEN
    CALL random_number(rand) ! [0-1]
    CALL quick_sort_asnd(rand,init,1,nbv)
  ELSE
    PRINT *, "  ERROR&SROP :: NOT SUCH OPTION in get_resampling_mtx for CC ", CC
    STOP
  END IF



  !==> generate resampling mxm matrix
IF( DG =='ON' )THEN
  inum(:)   = 0
  !(1) :: gen selected paarticle
  DO j=1,nbv ! jth column
    DO i=1,nbv-1
      IF( rand(j)<=acc(i) )THEN
        inum(j) = i
        GOTO 10
      ENDIF
    END DO
    inum(j) = nbv
10  CONTINUE
  END DO
  !!print *, inum(:)

  pmat(:,:) = 0.0d0
  DO i=1,nbv ! (1) diagonal component
    k = inum(i)
    IF( k/=0 .and. pmat(k,k)<0.5d0 )THEN
      pmat(k,k) = 1.0d0
      inum(i)   = 0
    END IF 
  END DO

  DO i=1,nbv ! (2) off-diagonal component
    k = inum(i)
    IF( k/=0 )THEN
      DO j=1,nbv
        IF( sum(pmat(:,j))<0.5 )THEN
          pmat(k,j) = 1.0d0
          GOTO 20
        END IF
      END DO
    END IF 
20  CONTINUE
  END DO
ELSE IF( DG =='OF' )THEN
  pmat(:,:) = 0.0d0
  DO j=1,nbv ! jth column
    DO i=1,nbv-1
      IF( rand(j)<=acc(i) )THEN
        pmat(i,j) = 1.0d0
        GOTO 30
      END IF
    END DO
    pmat(nbv,j) = 1.0d0
30  CONTINUE
  END DO
ENDIF
END SUBROUTINE get_resampling_mtx

!=========================================================
recursive subroutine quick_sort_asnd(var,init,first,last)
!=========================================================
implicit none
  integer :: first, last, i, j, it
  integer :: init(*)
  real(8) :: var(*) , x,t

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
