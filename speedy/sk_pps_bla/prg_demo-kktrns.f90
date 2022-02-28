program main
  implicit none

  !===> setting of speedy
  integer, parameter :: nbv    = 40
  integer, parameter :: ncase  = 1     ! the product (1000)
!  integer, parameter :: ntime  = 10000 ! the product (10000)
  integer, parameter :: ntime  = 200 ! the product (10000)

  integer :: init(nbv,nbv)
  real(8) :: rand(nbv), pfwgh(nbv), acc(nbv)
  real(8) :: pmat_answ(nbv,nbv), pmat_time(nbv,nbv,ntime), pmat_kkkk(nbv,nbv,ntime), tmp8(nbv,nbv)
  real(8) ::                                               pvec_kkkk(nbv,    ntime) 
  real(8) :: aerr(ncase,ntime) , perr(ncase,ntime)
  
  !===> global variables
  integer :: i, j, k, l, itime, icase
  real(8) :: p_ave, p_std, m_ave, m_std, timer, peff
  character(255) :: fname

  integer t0, t1, t_rate, t_max, diff
!-----------------------------------------------------------------------
  aerr(:,:) = 0.0d0

  call system_clock(t0)
  DO ICASE=1,NCASE
    call random_number(rand)
    pfwgh(:) = rand(:) / sum( rand(:) )
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

    pmat_time(:,:,:) = 0.0d0
    pmat_kkkk(:,:,:) = 0.0d0
    DO itime=1,ntime ! MONTE-Carlo Resampling M times
      if( mod(itime,ntime)==0 ) then
        call system_clock(t1, t_rate, t_max)
        timer = dble(t1-t0) / dble(t_rate)
        print '(i6.6,a,i6.6,a,f15.2)', icase," / ", ncase,   "  timer since start (s) :: ", timer
      end if
      CALL get_resampling_mtx(nbv,acc,pmat_time(1:nbv,1:nbv,itime))
      
      !SLOW!DO j=1,nbv ; DO i=1,nbv
      !SLOW!  pmat_kkkk(i,j,itime) = sum( pmat_time(i,j,1:itime) ) / dble( itime )
      !SLOW!END DO     ; END DO
      
      IF( itime==1 )THEN
        pmat_kkkk(1:nbv,1:nbv,itime) = pmat_time(1:nbv,1:nbv,itime)
      ELSE
        pmat_kkkk(1:nbv,1:nbv,itime) = ( pmat_kkkk(1:nbv,1:nbv,itime-1)*(dble(itime)-1.0d0) + pmat_time(1:nbv,1:nbv,itime) ) / dble( itime )
      END IF
      DO i=1,nbv
        pvec_kkkk(i,itime)           = sum( pmat_kkkk(i,1:nbv,itime) ) / dble( nbv )
      END DO
      !!DO i=1,nbv
      !!  print '(40i1)', int( pmat_time(i,1:nbv,itime) )
      !!END DO
    END DO

    DO j=1,nbv ; DO i=1,nbv
      pmat_answ(i,j) = sum( pmat_time(i,j,1:ntime) ) / dble( ntime )
    END DO     ; END DO

    DO itime=1,ntime
      aerr(icase,itime) = sum(  abs(pmat_answ(:,:)-pmat_kkkk(:,:,itime)) )
      perr(icase,itime) = sum(  abs(pfwgh(:)      -pvec_kkkk(:  ,itime)) )
    END DO

    IF( icase==1 )THEN
      write(fname,'(a,i6.6,a,i8.8,a,i8.8,a)')"./out_900/sampled-error_M",nbv,"NTIME",ntime,"NCASE",ncase,".bin"
      open(1,file=trim(fname),form="unformatted",access="direct",recl=nbv*nbv*4,action="write")
        write(1,rec=1) (( sngl( tmp8(i,j)            ),j=1,nbv),i=1,nbv)
        write(1,rec=2) (( sngl( pmat_time(i,j,1)     ),j=1,nbv),i=1,nbv)
        write(1,rec=3) (( sngl( pmat_time(i,j,2)     ),j=1,nbv),i=1,nbv)
        write(1,rec=4) (( sngl( pmat_kkkk(i,j,200)   ),j=1,nbv),i=1,nbv)
        write(1,rec=5) (( sngl( pmat_kkkk(i,j,ntime) ),j=1,nbv),i=1,nbv)
      close(1)
    ENDIF
  END DO

    print '(a)', "========================ave time ============================="
    do i=1,nbv
      print '(I3,a,3f8.3,a,3f8.3)', i, &
        "   || one-su, mlt-su, ans :: ",sum(pmat_time(i,:,1)), sum(pmat_kkkk(i,:,200)), sum(pmat_kkkk(i,:,ntime)), &
        "   || wgh, acc, peff :: ",pfwgh(i)*dble(nbv), acc(i), peff
    end do

  !==> statistics
  write(fname,'(a,i6.6,a,i8.8,a,i8.8,a)')"./out_900/sampled-error_M",nbv,"NTIME",ntime,"NCASE",ncase,".txt"
  open(1,file=trim(fname),form="formatted",action="write")
  DO itime=1,ntime
    p_ave = sum( perr(1:ncase,itime) ) / dble( NCASE )
    p_std = dsqrt( sum( (perr(1:ncase,itime)-p_ave)**2.0d0 ) / dble( NCASE ) )
    m_ave = sum( aerr(1:ncase,itime) ) / dble( NCASE )
    m_std = dsqrt( sum( (aerr(1:ncase,itime)-m_ave)**2.0d0 ) / dble( NCASE ) )

    write(1,'(2i,6f15.8,4x,6f15.8)') ncase, itime,                                                          & 
      p_ave, p_std, p_ave-p_std, p_ave+p_std, minval( perr(1:ncase,itime) ), maxval( perr(1:ncase,itime) ), &
      m_ave, m_std, m_ave-m_std, m_ave+m_std, minval( aerr(1:ncase,itime) ), maxval( aerr(1:ncase,itime) )
    !!write(6,'(2i,4f15.8)') ncase, itime, ave, std, ave-std, ave+std
  END DO
  close(1)

  print '(a)', "program end from prg_demo-kktrns.f90"
  stop
end

!-----------------------------------------------------------------------
SUBROUTINE get_resampling_mtx(nbv,acc,pmat)
!-----------------------------------------------------------------------
IMPLICIT NONE
  INTEGER, INTENT(in)  :: nbv
  REAL(8), INTENT(in)  :: acc(nbv)
  REAL(8), INTENT(out) :: pmat(nbv,nbv)

  INTEGER :: i, j, k, init(nbv), inum(nbv)
  REAL(8) :: rand(nbv)

  DO j=1,nbv ; init(j) = j  ; END DO
  CALL random_number(rand) ! [0-1]
  CALL quick_sort_asnd(rand,init,1,nbv)

  !==> generate resampling mxm matrix
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

!ORG!  pmat(:,:) = 0.0d0
!ORG!  DO j=1,nbv ! jth column
!ORG!    DO i=1,nbv-1
!ORG!      IF( rand(j)<=acc(i) )THEN
!ORG!        pmat(i,j) = 1.0d0
!ORG!        GOTO 10
!ORG!      END IF
!ORG!    END DO
!ORG!    pmat(nbv,j) = 1.0d0
!ORG!10  CONTINUE
!ORG!  END DO
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
