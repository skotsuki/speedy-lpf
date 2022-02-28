MODULE filter
!=======================================================================
!
! [PURPOSE:] Module for LETKF with SPEEDY
!
! [HISTORY:]
!   06/XX 2018 P. Andrew  :: /data6/apensoneault/TEST-LETKF-EXPERIMENT/speedy/backup_letkf_7_26/old/filter5.f90
!   04/21 2019 S. Kotsuki :: edits for SPEEDY-LPF
!
!=======================================================================
  USE common
  USE common_mpi
  USE common_speedy
  USE common_mpi_speedy
  USE common_letkf
  USE letkf_obs
  IMPLICIT NONE

  PRIVATE
  !SK201904!PUBLIC :: Fr_mat, Fr_vec,  LANCZOS_FILTER_mat, LANCZOS_FILTER_vec,LANCZOS_FILTER_mat_naive, LANCZOS_FILTER_vec_naive
  PUBLIC :: Fr_mat, Fr_vec,  LANCZOS_FILTER_mat, LANCZOS_FILTER_vec

  INTEGER, PARAMETER :: ISC=1
  INTEGER, PARAMETER :: NTRUN=30, MTRUN=30, IX=96, IY=24 
  INTEGER, PARAMETER :: NX=NTRUN+2, MX=MTRUN+1 , MXNX=MX*NX,MX2=2*MX
  INTEGER, PARAMETER :: IL=2*IY, NTRUN1=NTRUN+1 
  INTEGER, PARAMETER :: NXP=NX+1 , MXP=ISC*MTRUN+1, LMAX=MXP+NX-2 
  INTEGER, parameter      :: a=20 !<-- AP 6/18/2018 - maximum number of coefficients
  REAL(r_size), parameter :: fx=1.0/5.0 !<-- AP 6/18/2018 - cuttoff frequency x (what should this be?)
  REAL(r_size), parameter :: fy=1.0/5.0 !<-- AP 6/19/2018 - cuttoff frequency y
  INTEGER, PARAMETER :: L=21, M=21, LDIM=nlon
CONTAINS
  
  REAL(r_size) FUNCTION TSINC(X,k)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: X
  INTEGER, INTENT(IN) :: k
  IF(X==0) THEN
    TSINC=1
  ELSEIF(k<a) then
    TSINC=SIN(X)/X
  ELSE
    TSINC=0
  END IF
  RETURN
  END FUNCTION TSINC

  REAL(r_size) FUNCTION SINC(X)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: X
  IF(X==0) THEN
    SINC=1
  ELSE
    SINC=SIN(X)/X
  END IF
  RETURN
  END FUNCTION SINC

 
  REAL(r_size) FUNCTION BESINC(X)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: X
  IF(X==0) THEN
    BESINC=1/2
  ELSE
    BESINC=BESSEL_JN(1,X)/X
  END IF
  RETURN
  END FUNCTION BESINC

  SUBROUTINE Fr_mat(trans,trans_filt)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: trans(nlon,nlat,nbv,nbv)
  REAL(r_size), INTENT(OUT) :: trans_filt(nlon,nlat,nbv,nbv)
  INTEGER :: i,j,IER,LENSAV,LENWRK
  real(r_size), allocatable, dimension ( : ) :: WORK
  real(r_size), allocatable, dimension ( : ) :: WSAVE
  COMPLEX(r_size) :: fourier_trans(nlon,nlat,nbv,nbv)
  LENSAV=MAX(2*(L+M) + INT(LOG(REAL(L))/LOG(2.)) + INT(LOG(REAL(M))/LOG(2.)) + 8,&
         L+int(log(real(L,kind = 8 ))/log(2.0D+00))+4+2*M+int(log(real(M,kind=8))/log(2.0D+00))+4+M+int(log(real(M,kind=8))/log(2.0D+00))+4)
  LENWRK=max(LDIM*M,M*(L+1),2*M*L)
  ALLOCATE(WSAVE(1:LENSAV),WORK(1:LENWRK))
  CALL cfft2i (L,M,WSAVE,LENSAV,IER)
  trans_filt=0.0d0
  fourier_trans=0.0d0
  do i=1,nbv
    do j=1,nbv
      fourier_trans(:,:,i,j)=cmplx(trans(:,:,i,j))
      CALL CFFT2F (LDIM, L, M, fourier_trans(:,:,i,j), WSAVE, LENSAV, WORK, LENWRK, IER) 
      CALL CFFT2B (LDIM, L, M, fourier_trans(:,:,i,j), WSAVE, LENSAV, WORK, LENWRK, IER)
      trans_filt(:,:,i,j)=REAL(fourier_trans(:,:,i,j),r_size)
    end do
  end do
  deallocate (WORK)
  deallocate (WSAVE)
  RETURN
  END SUBROUTINE Fr_mat

  SUBROUTINE Fr_vec(vec,vec_filt)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: vec(nlon,nlat,nbv)
  REAL(r_size), INTENT(OUT) :: vec_filt(nlon,nlat,nbv)
  INTEGER :: i,IER,LENSAV,LENWRK
  real(r_size), allocatable, dimension ( : ) :: WORK
  real(r_size), allocatable, dimension ( : ) :: WSAVE
  COMPLEX(r_size) :: fourier_vec(nlon,nlat,nbv)
  LENSAV=MAX(2*(L+M) + INT(LOG(REAL(L))/LOG(2.)) + INT(LOG(REAL(M))/LOG(2.)) + 8,&
         L+int(log(real(L,kind = 8 ))/log(2.0D+00))+4+2*M+int(log(real(M,kind=8))/log(2.0D+00))+4+M+int(log(real(M,kind=8))/log(2.0D+00))+4)
  LENWRK=max(LDIM*M,M*(L+1),2*M*L)
  ALLOCATE(WSAVE(1:LENSAV),WORK(1:LENWRK))
  CALL cfft2i (L,M,WSAVE,LENSAV,IER)
  vec_filt=0.0d0
  fourier_vec=0.0d0
  do i=1,nbv
    fourier_vec(:,:,i)=cmplx(vec(:,:,i))
    CALL CFFT2F (LDIM, L, M, fourier_vec(:,:,i), WSAVE, LENSAV, WORK, LENWRK, IER)
    CALL CFFT2B (LDIM, L, M, fourier_vec(:,:,i), WSAVE, LENSAV, WORK, LENWRK, IER)
    vec_filt(:,:,i)=REAL(fourier_vec(:,:,i),r_size)
  end do
  deallocate (WORK)
  deallocate (WSAVE)
  RETURN
  END SUBROUTINE Fr_vec





!=======================================================================================
  SUBROUTINE LANCZOS_FILTER_mat(trans,trans_filt)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: trans(nlon,nlat,nbv,nbv)
  REAL(r_size), INTENT(OUT) :: trans_filt(nlon,nlat,nbv,nbv)
  REAL(r_size) :: u,upp,ump,umm,upm
  INTEGER :: i,j,k,n,o
  INTEGER :: IER,LENSAV,LENWRK
  real(r_size), allocatable, dimension ( : ) :: WORK
  real(r_size), allocatable, dimension ( : ) :: WSAVE
  COMPLEX(r_size) :: fourier_trans(nlon,nlat,nbv,nbv)
  COMPLEX(r_size) :: fourier_filt(nlon,nlat,nbv,nbv)
  real(r_size) :: weight(nlon,nlat)
  COMPLEX(r_size) :: fourier_weight(nlon,nlat)
  LENSAV=MAX(2*(L+M) + INT(LOG(REAL(L))/LOG(2.)) + INT(LOG(REAL(M))/LOG(2.)) + 8,&
         L+int(log(real(L,kind = 8 ))/log(2.0D+00))+4+2*M+int(log(real(M,kind=8))/log(2.0D+00))+4+M+int(log(real(M,kind=8))/log(2.0D+00))+4)
  LENWRK=max(LDIM*M,M*(L+1),2*M*L)
  ALLOCATE(WSAVE(1:LENSAV),WORK(1:LENWRK))
  CALL cfft2i (L,M,WSAVE,LENSAV,IER)
  trans_filt=0.0d0
  fourier_filt=0.0d0
  fourier_trans=0.0d0

  do i=1,nbv
    do j=1,nbv
      fourier_trans(:,:,i,j)=cmplx(trans(:,:,i,j))
      CALL CFFT2F (LDIM, L, M, fourier_trans(:,:,i,j), WSAVE, LENSAV, WORK, LENWRK, IER)
    end do
  end do

!  do k=0,M/2
!   do n=0,L/2
!     uplus=sqrt((fx*n)**2+(fy*k)**2)
!     uminus
!     weight(n+1,k+1)=real(L,r_size)*real(M,r_size)*fx*fy*SINC(PI*real(k)/real(M))*SINC(PI*REAL(n)/REAL(L))*BESSEL_JN(1,2*PI*u)/u
!     if(k<M/2) then
!       if (k<M/2) then
!        weight(n+1,k+1)=4*M*L*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M/2))*sinc(pi*n/(L/2))
!       end if
!     end if
!     weight(L-n,M-k)=4*M*L*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M/2))*sinc(pi*n/(L/2))
!   end do
! end do  
!  do n=0,M-1
!    do k=0,L-1
!      u=sqrt((fx*n)**2+(fy*k)**2)
!      weight(n+1,k+1)=fx*fy*SINC(PI*real(k)/real(M))*SINC(PI*REAL(n)/REAL(L))*2*PI*BESINC(2*PI*u)
!    end do
!  end do

  do n=0,M-1
    do k=0,L-1
      upp=sqrt((fx*(n))**2+(fy*(k))**2)
      weight(n+1,k+1)=fx*fy*TSINC(PI*real(k)/real(M),k)*TSINC(PI*REAL(n)/REAL(L),n)*2*PI*BESINC(2*PI*upp)
    end do
  end do
  fourier_weight(:,:)=cmplx(weight(:,:))
  CALL CFFT2F (LDIM, L, M,fourier_weight, WSAVE, LENSAV, WORK, LENWRK, IER)


!  do n=0,(M/2-1)
!    do k=0,(L/2-1)
!      upp=sqrt((fx*n)**2+(fy*k)**2)
!      upm=sqrt((fx*n)**2+(fy*(k-M/2))**2)
!      ump=sqrt((fx*(n-L/2))**2+(fy*k)**2)
!      umm=sqrt((fx*(n-L/2))**2+(fy*(k-M/2))**2)
!      weight(n+1,k+1)=fx*fy*SINC(PI*k/M)*SINC(PI*n/L)*2*PI*BESINC(2*PI*upp)
!      weight(n+1,M/2-k)=fx*fy*SINC(PI*(k-M/2)/M)*SINC(PI*n/L)*2*PI*BESINC(2*PI*upm)
!      weight(L/2-n,k+1)=fx*fy*SINC(PI*k/M)*SINC(PI*(n-L/2)/L)*2*PI*BESINC(2*PI*ump)
!      weight(L/2-n,M/2-k)=fx*fy*SINC(PI*(k-M/2)/M)*SINC(PI*(n-L/2)/L)*2*PI*BESINC(2*PI*umm)
!    end do
!  end do
 
  do i=1,nbv
    do j=1,nbv
      fourier_filt(:,:,i,j)=fourier_trans(:,:,i,j)*fourier_weight(:,:)/abs(fourier_weight(1,1))
      CALL CFFT2B (LDIM, L, M, fourier_filt(:,:,i,j), WSAVE, LENSAV, WORK, LENWRK, IER)
      trans_filt(:,:,i,j)=REAL(fourier_filt(:,:,i,j),r_size)
    end do
  end do
  !PRINT *, maxval(abs(trans_filt(:,:,:,:)-trans(:,:,:,:))), minval(abs(trans_filt(:,:,:,:)-trans(:,:,:,:)))


  deallocate (WORK)
  deallocate (WSAVE)
  RETURN
  END SUBROUTINE LANCZOS_FILTER_mat

!=======================================================================================
  SUBROUTINE LANCZOS_FILTER_vec(vec,vec_filt)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: vec(nlon,nlat,nbv)
  REAL(r_size), INTENT(OUT) :: vec_filt(nlon,nlat,nbv)
  COMPLEX(r_size) :: fourier_vec(nlon,nlat,nbv)
  COMPLEX(r_size) :: fourier_filt(nlon,nlat,nbv)
  REAL(r_size) :: u,upp,ump,umm,upm
  INTEGER :: i,j,k,n,o
  INTEGER :: IER,LENSAV,LENWRK
  real(r_size), allocatable, dimension ( : ) :: WORK
  real(r_size), allocatable, dimension ( : ) :: WSAVE
  real(r_size) :: weight(nlon,nlat)
  COMPLEX(r_size) :: fourier_weight(nlon,nlat)
  LENSAV=max(2*(L+M) + INT(LOG(REAL(L))/LOG(2.)) + INT(LOG(REAL(M))/LOG(2.)) + 8,L+int(log(real(L,kind = 8 ))/log(2.0D+00))+4+2*M&
             +int(log(real(M,kind=8))/log(2.0D+00))+4+M+int(log(real(M,kind=8))/log(2.0D+00))+4)
  LENWRK=max(LDIM*M,M*(L+1),2*L*M)
  ALLOCATE(WSAVE(1:LENSAV),WORK(1:LENWRK))
  CALL cfft2i (L,M,WSAVE,LENSAV,IER)
  vec_filt=0.0d0
  fourier_vec=cmplx(0.0d0)
  fourier_filt=cmplx(0.0d0)

  do i=1,nbv
    fourier_vec(:,:,i)=cmplx(vec(:,:,i))
    CALL CFFT2F (LDIM, L, M,fourier_vec(:,:,i) , WSAVE, LENSAV, WORK, LENWRK, IER)
  end do

!  do n=1,M/2
!    do k=1,L/2
!      upp=sqrt((fx*(n-1))**2+(fy*(k-1))**2)
!      upm=sqrt((fx*(n-1))**2+(fy*(M/2-k))**2)
!      ump=sqrt((fx*(L/2-n))**2+(fy*(k-1))**2)
!      umm=sqrt((fx*(L/2-n))**2+(fy*(M/2-k))**2)
!      weight(n,k)=fx*fy*SINC(PI*real(k-1)/real(M))*SINC(PI*REAL(n-1)/REAL(L))*2*PI*BESINC(2*PI*upp)
!      weight(n,M/2-k)=fx*fy*SINC(PI*real(M/2-k)/real(M))*SINC(PI*REAL(n-1)/REAL(L))*2*PI*BESINC(2*PI*upm)
!      weight(L/2-n,k)=fx*fy*SINC(PI*real(k-1)/real(M))*SINC(PI*REAL(L/2-n)/REAL(L))*2*PI*BESINC(2*PI*ump)
!      weight(L/2-n,M/2-k)=fx*fy*SINC(PI*real(M/2-k)/real(M))*SINC(PI*REAL(L/2-n)/REAL(L))*2*PI*BESINC(2*PI*umm)
!    end do
!  end do

!!!  do n=1,(M/2-1)
!!!    do k=1,(L/2-1)
!!!      upp=sqrt((fx*(n-1))**2+(fy*(k-1))**2)
!!!      upm=sqrt((fx*(n-1))**2+(fy*(M/2-k))**2)
!!!      ump=sqrt((fx*(L/2-n))**2+(fy*(k-1))**2)
!!!      umm=sqrt((fx*(L/2-n))**2+(fy*(M/2-k))**2)
!!!      weight(n,k)=fx*fy*SINC(PI*(k-1)/M)*SINC(PI*(n-1)/L)*2*PI*BESINC(2*PI*upp)
!!!      weight(n,M/2-k)=fx*fy*SINC(PI*(M/2-k)/M)*SINC(PI*(n-1)/L)*2*PI*BESINC(2*PI*upm)
!!!      weight(L/2-n,k)=fx*fy*SINC(PI*(k-1)/M)*SINC(PI*(L/2-n)/L)*2*PI*BESINC(2*PI*ump)
!!!      weight(L/2-n,M/2-k)=fx*fy*SINC(PI*(M/2-k)/M)*SINC(PI*(L/2-n)/L)*2*PI*BESINC(2*PI*umm)
!!!    end do
!!!  end do

!  do n=0,(M/2-1)
!    do k=0,(L/2-1)
!      upp=sqrt((fx*n)**2+(fy*k)**2)
!      upm=sqrt((fx*n)**2+(fy*(k-M/2))**2)
!      ump=sqrt((fx*(n-L/2))**2+(fy*k)**2)
!      umm=sqrt((fx*(n-L/2))**2+(fy*(k-M/2))**2)
!      weight(n+1,k+1)=fx*fy*SINC(PI*k/M)*SINC(PI*n/L)*2*PI*BESINC(2*PI*upp)
!      weight(n+1,M/2-k)=fx*fy*SINC(PI*(k-M/2)/M)*SINC(PI*n/L)*2*PI*BESINC(2*PI*upm)
!      weight(L/2-n,k+1)=fx*fy*SINC(PI*k/M)*SINC(PI*(n-L/2)/L)*2*PI*BESINC(2*PI*ump)
!      weight(L/2-n,M/2-k)=fx*fy*SINC(PI*(k-M/2)/M)*SINC(PI*(n-L/2)/L)*2*PI*BESINC(2*PI*umm)
!    end do
!  end do


  do n=0,M-1
    do k=0,L-1
      upp=sqrt((fx*(n))**2+(fy*(k))**2)
      weight(n+1,k+1)=fx*fy*TSINC(PI*real(k)/real(M),k)*TSINC(PI*REAL(n)/REAL(L),n)*2*PI*BESINC(2*PI*upp)
    end do
  end do
  


  fourier_weight(:,:)=cmplx(weight(:,:))

  CALL CFFT2F (LDIM, L, M,fourier_weight, WSAVE, LENSAV, WORK, LENWRK, IER)
  do i=1,nbv
    fourier_filt(:,:,i)=fourier_vec(:,:,i)*fourier_weight(:,:)/abs(fourier_weight(1,1))
    CALL CFFT2B (LDIM, L, M, fourier_filt(:,:,i), WSAVE, LENSAV, WORK, LENWRK, IER)
    vec_filt(:,:,i)=REAL(fourier_filt(:,:,i),r_size)
  end do
  deallocate (WORK)
  deallocate (WSAVE)
  RETURN
  END SUBROUTINE LANCZOS_FILTER_vec


!SK201904!  SUBROUTINE LANCZOS_FILTER_mat_naive(trans,trans_filt)
!SK201904!  IMPLICIT NONE
!SK201904!  REAL(r_size), INTENT(IN) :: trans(nlon,nlat,nbv,nbv)
!SK201904!  REAL(r_size), INTENT(OUT) :: trans_filt(nlon,nlat,nbv,nbv)
!SK201904!  REAL(r_size) :: u
!SK201904!  INTEGER :: i,j,k,n,o
!SK201904!  INTEGER :: IER,LENSAV,LENWRK
!SK201904!  real(r_size), allocatable, dimension ( : ) :: WORK
!SK201904!  real(r_size), allocatable, dimension ( : ) :: WSAVE
!SK201904!  COMPLEX(r_size) :: fourier_trans(nlon,nlat,nbv,nbv)
!SK201904!  COMPLEX(r_size) :: fourier_filt(nlon,nlat,nbv,nbv)
!SK201904!  real(r_size) :: weight(nlon,nlat)
!SK201904!  COMPLEX(r_size) :: fourier_weight(nlon,nlat)
!SK201904!  LENSAV=MAX(2*(L+M) + INT(LOG(REAL(L))/LOG(2.)) + INT(LOG(REAL(M))/LOG(2.)) + 8,&
!SK201904!         L+int(log(real(L,kind = 8 ))/log(2.0D+00))+4+2*M+int(log(real(M,kind=8))/log(2.0D+00))+4+M+int(log(real(M,kind=8))/log(2.0D+00))+4)
!SK201904!  LENWRK=max(LDIM*M,M*(L+1),2*M*L)
!SK201904!  ALLOCATE(WSAVE(1:LENSAV),WORK(1:LENWRK)) 
!SK201904!  CALL cfft2i (L,M,WSAVE,LENSAV,IER)
!SK201904!  trans_filt=0.0d0
!SK201904!  fourier_filt=0.0d0
!SK201904!  fourier_trans=0.0d0
!SK201904!  
!SK201904!  do i=1,nbv
!SK201904!    do j=1,nbv
!SK201904!      fourier_trans(:,:,i,j)=cmplx(trans(:,:,i,j))
!SK201904!      CALL CFFT2F (LDIM, L, M, fourier_trans(:,:,i,j), WSAVE, LENSAV, WORK, LENWRK, IER)
!SK201904!    end do
!SK201904!  end do
!SK201904!!  do k=0,M/2
!SK201904!!    do n=0,L/2
!SK201904!!      if(k<M/2) then
!SK201904!!        if (k<M/2) then
!SK201904!!          weight(n+1,k+1)=4*M*L*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M/2))*sinc(pi*n/(L/2))
!SK201904!!        end if
!SK201904!!      end if
!SK201904!!      weight(L-n,M-k)=4*M*L*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M/2))*sinc(pi*n/(L/2))
!SK201904!!    end do
!SK201904!!  end do
!SK201904!! do k=1,M/2
!SK201904!!   do n=1,L/2
!SK201904!!    weight(n,k)=4*M*L*fy*fx*sinc(2*PI*fx*(k-1))*sinc(2*PI*fy*(n-1))*sinc(pi*(k-1)/(M/2))*sinc(pi*(n-1)/(L/2))
!SK201904!!    weight(L-n,k)=4*M*L*fy*fx*sinc(2*PI*fx*(k-1))*sinc(2*PI*fy*(-n))*sinc(pi*(k-1)/(M/2))*sinc(pi*(-n)/(L/2))
!SK201904!!    weight(n,M-k)=4*M*L*fy*fx*sinc(2*PI*fx*(-k))*sinc(2*PI*fy*(n-1))*sinc(pi*(-k)/(M/2))*sinc(pi*(n-1)/(L/2))
!SK201904!!    weight(L-n,M-k)=4*M*L*fy*fx*sinc(2*PI*fx*(-k))*sinc(2*PI*fy*(-n))*sinc(pi*(-k)/(M/2))*sinc(pi*(-n)/(L/2))
!SK201904!!    weight(n+1,k+1)=4*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M))*sinc(pi*n/(L))
!SK201904!!   end do
!SK201904!! end do
!SK201904!
!SK201904!!  do k=0,M/2
!SK201904!!    do n=0,L/2
!SK201904!!      if(k<M/2) then
!SK201904!!        if (n<L/2) then
!SK201904!!          weight(n+1,k+1)=4*M*L*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M/2))*sinc(pi*n/(L/2))
!SK201904!!        end if
!SK201904!!        weight(n+1,L-k)=4*M*L*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M/2))*sinc(pi*n/(L/2))
!SK201904!!      end if
!SK201904!!      if (n<L/2) then
!SK201904!!        weight(L-n,k+1)=4*M*L*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M/2))*sinc(pi*n/(L/2))
!SK201904!!      end if
!SK201904!!      weight(L-n,M-k)=4*M*L*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M/2))*sinc(pi*n/(L/2))
!SK201904!!    end do
!SK201904!!  end do
!SK201904!
!SK201904!
!SK201904!! do k=0,M-1
!SK201904!!    do n=0,L-1
!SK201904!!      u=sqrt((fx*n)**2+(fy*k)**2)
!SK201904!!      !weight(n+1,k+1)=real(L,r_size)*real(M,r_size)*fx*fy*SINC(PI*real(k)/real(M))*SINC(PI*REAL(n)/REAL(L))*BESSEL_JN(1,2*PI*u)/u
!SK201904!!      weights(n+1,k+1)=4*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/M)*sinc(pi*n/L)
!SK201904!!    end do
!SK201904!! end do
!SK201904!  !weight(1,1)=fx*fy*PI
!SK201904! 
!SK201904!  do k=0,M-1
!SK201904!    do n=0,L-1
!SK201904!      weight(n+1,k+1)=4*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/M)*sinc(pi*n/L)
!SK201904!    end do
!SK201904!  end do
!SK201904!  fourier_weight(:,:)=cmplx(weight(:,:))
!SK201904!  CALL CFFT2F (LDIM, L, M,fourier_weight, WSAVE, LENSAV, WORK, LENWRK, IER)
!SK201904!
!SK201904!  do i=1,nbv 
!SK201904!    do j=1,nbv
!SK201904!      fourier_filt(:,:,i,j)=fourier_trans(:,:,i,j)*fourier_weight(:,:)
!SK201904!      CALL CFFT2B (LDIM, L, M, fourier_filt(:,:,i,j), WSAVE, LENSAV, WORK, LENWRK, IER)
!SK201904!      trans_filt(:,:,i,j)=REAL(fourier_filt(:,:,i,j),r_size) 
!SK201904!    end do
!SK201904!  end do
!SK201904!
!SK201904!!  do i=1,L
!SK201904!!    do j=1,M
!SK201904!!      do k=-a,a
!SK201904!!        do n=-a,-1
!SK201904!!          u=sqrt((fx*n)**2+(fy*k)**2)
!SK201904!!          fourier_filt(i,j,:,:)=fourier_filt(i,j,:,:)+exp(cmplx(0.,1.)*2*PI*((REAL((i-1)*n)/REAL(nlon-1))+(REAL((j-1)*k)/REAL(nlat-1))))&
!SK201904!!                              *fx*fy*SINC(PI*real(k)/real(a))*SINC(PI*REAL(n)/REAL(a))*BESSEL_JN(1,2*PI*u)/u
!SK201904!!        end do
!SK201904!!        do n=1,a
!SK201904!!          u=sqrt((fx*n)**2+(fy*k)**2)
!SK201904!!          fourier_filt(i,j,:,:)=fourier_filt(i,j,:,:)+exp(cmplx(0.,1.)*2*PI*((REAL((i-1)*n)/REAL(nlon-1))+(REAL((j-1)*k)/REAL(nlat-1))))&
!SK201904!!                              *fx*fy*SINC(PI*real(k)/real(a))*SINC(PI*REAL(n)/REAL(a))*BESSEL_JN(1,2*PI*u)/u
!SK201904!!        end do
!SK201904!!        n=0
!SK201904!!        if(k==0) then
!SK201904!!          u=PI
!SK201904!!          fourier_filt(i,j,:,:)=fourier_filt(i,j,:,:)+fx*fy*PI
!SK201904!!        else
!SK201904!!          u=sqrt((fx*n)**2+(fy*k)**2)
!SK201904!!          u= BESSEL_JN(1,2*PI*u)/u
!SK201904!!          fourier_filt(i,j,:,:)=fourier_filt(i,j,:,:)+exp(cmplx(0.,1.)*2*PI*(REAL((j-1)*k)/REAL(nlat-1)))*fx*fy*SINC(PI*real(k)/real(a))*u
!SK201904!!        end if
!SK201904!!      end do
!SK201904!!      fourier_filt(i,j,:,:)=fourier_filt(i,j,:,:)*fourier_trans(i,j,:,:)
!SK201904!!    end do
!SK201904!!  end do
!SK201904!!  do i=1,nbv
!SK201904!!    do j=1,nbv
!SK201904!!      CALL CFFT2B (LDIM, L, M, fourier_filt(:,:,i,j), WSAVE, LENSAV, WORK, LENWRK, IER)
!SK201904!!      trans_filt(:,:,i,j)=REAL(fourier_filt(:,:,i,j),r_size)
!SK201904!!    end do
!SK201904!!  end do
!SK201904!  deallocate (WORK)
!SK201904!  deallocate (WSAVE)
!SK201904!  RETURN
!SK201904!  END SUBROUTINE LANCZOS_FILTER_mat_naive
!SK201904!
!SK201904!
!SK201904!  SUBROUTINE LANCZOS_FILTER_vec_naive(vec,vec_filt)
!SK201904!  IMPLICIT NONE
!SK201904!  REAL(r_size), INTENT(IN) :: vec(nlon,nlat,nbv)
!SK201904!  REAL(r_size), INTENT(OUT) :: vec_filt(nlon,nlat,nbv)
!SK201904!  COMPLEX(r_size) :: fourier_vec(nlon,nlat,nbv)
!SK201904!  COMPLEX(r_size) :: fourier_filt(nlon,nlat,nbv)
!SK201904!  REAL(r_size) :: u
!SK201904!  INTEGER :: i,j,k,n,o
!SK201904!  INTEGER :: IER,LENSAV,LENWRK
!SK201904!  real(r_size), allocatable, dimension ( : ) :: WORK
!SK201904!  real(r_size), allocatable, dimension ( : ) :: WSAVE
!SK201904!  real(r_size) :: weight(nlon,nlat)
!SK201904!  COMPLEX(r_size) :: fourier_weight(nlon,nlat)
!SK201904!  LENSAV=max(2*(L+M) + INT(LOG(REAL(L))/LOG(2.)) + INT(LOG(REAL(M))/LOG(2.)) + 8,L+int(log(real(L,kind = 8 ))/log(2.0D+00))+4+2*M&
!SK201904!             +int(log(real(M,kind=8))/log(2.0D+00))+4+M+int(log(real(M,kind=8))/log(2.0D+00))+4)
!SK201904!  LENWRK=max(LDIM*M,M*(L+1),2*L*M)
!SK201904!  ALLOCATE(WSAVE(1:LENSAV),WORK(1:LENWRK))
!SK201904!  CALL cfft2i (L,M,WSAVE,LENSAV,IER)
!SK201904!  vec_filt=0.0d0
!SK201904!  fourier_vec=cmplx(0.0d0)
!SK201904!  fourier_filt=cmplx(0.0d0)
!SK201904!
!SK201904!  do i=1,nbv
!SK201904!    fourier_vec(:,:,i)=cmplx(vec(:,:,i))
!SK201904!    CALL CFFT2F (LDIM, L, M,fourier_vec(:,:,i) , WSAVE, LENSAV, WORK, LENWRK, IER) 
!SK201904!  end do
!SK201904!
!SK201904!  do k=0,M-1
!SK201904!    do n=0,L-1
!SK201904!      weight(n+1,k+1)=4*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M))*sinc(pi*n/(L))
!SK201904!    end do
!SK201904!  end do
!SK201904!
!SK201904!
!SK201904!
!SK201904! 
!SK201904!! do k=1,M/2
!SK201904!!   do n=1,L/2
!SK201904!!    ! if(k<M/2) then
!SK201904!!    !   if (n<L/2) then
!SK201904!!    !     weight(n+1,k+1)=4*M*L*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M/2))*sinc(pi*n/(L/2))
!SK201904!!    !   end if
!SK201904!!    !   weight(L-n,k+1)=4*M*L*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M/2))*sinc(pi*n/(L/2))
!SK201904!!    ! end if
!SK201904!!    ! if (n<L/2) then
!SK201904!!    !   weight(L-n,k+1)=4*M*L*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M/2))*sinc(pi*n/(L/2))
!SK201904!!    ! end if
!SK201904!!    ! weight(L-n,M-k)=4*M*L*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M/2))*sinc(pi*n/(L/2))
!SK201904!!    !weight(n,k)=4*M*L*fy*fx*sinc(2*PI*fx*(k-1))*sinc(2*PI*fy*(n-1))*sinc(pi*(k-1)/(M/2))*sinc(pi*(n-1)/(L/2))
!SK201904!!    !weight(L-n,k)=4*M*L*fy*fx*sinc(2*PI*fx*(k-1))*sinc(2*PI*fy*(-n))*sinc(pi*(k-1)/(M/2))*sinc(pi*(-n)/(L/2))
!SK201904!!    !weight(n,M-k)=4*M*L*fy*fx*sinc(2*PI*fx*(-k))*sinc(2*PI*fy*(n-1))*sinc(pi*(-k)/(M/2))*sinc(pi*(n-1)/(L/2))
!SK201904!!    !weight(L-n,M-k)=4*M*L*fy*fx*sinc(2*PI*fx*(-k))*sinc(2*PI*fy*(-n))*sinc(pi*(-k)/(M/2))*sinc(pi*(-n)/(L/2))
!SK201904!!    !weight(n+1,k+1)=4*fy*fx*sinc(2*PI*fx*k)*sinc(2*PI*fy*n)*sinc(pi*k/(M))*sinc(pi*n/(L))
!SK201904!!   end do
!SK201904!! end do
!SK201904!  !weight(1,1)=fx*fy*PI
!SK201904!  fourier_weight(:,:)=cmplx(weight(:,:))
!SK201904! 
!SK201904!  CALL CFFT2F (LDIM, L, M,fourier_weight, WSAVE, LENSAV, WORK, LENWRK, IER)
!SK201904!  do i=1,nbv
!SK201904!    fourier_filt(:,:,i)=fourier_vec(:,:,i)*fourier_weight(:,:)
!SK201904!    CALL CFFT2B (LDIM, L, M, fourier_filt(:,:,i), WSAVE, LENSAV, WORK, LENWRK, IER)
!SK201904!    vec_filt(:,:,i)=REAL(fourier_filt(:,:,i),r_size)
!SK201904!  end do
!SK201904!
!SK201904!!  PRINT *, 'UNFILTERED+++++++++++++++++++++++++++++++++++++++++++++++'
!SK201904!!  PRINT *, vec(:,:,1)-vec_filt(:,:,1)
!SK201904!!  PRINT *, 'FILTERED*************************************************'
!SK201904!!  PRINT *, 'UNFILTERED+++++++++++++++++++++++++++++++++++++++++++++++'
!SK201904!!  PRINT *, vec(:,:,1)
!SK201904!!  PRINT *, 'FILTERED*************************************************'
!SK201904!!  PRINT *, fourier_filt(:,:,1)
!SK201904!!  PRINT *, '_________________________________________________________'
!SK201904!!  do i=1,L
!SK201904!!    do j=1,M
!SK201904!!      do k=-a,a
!SK201904!!        do n=-a,-1
!SK201904!!          u=sqrt((fx*n)**2+(fy*k)**2)
!SK201904!!          fourier_filt(i,j,:)=fourier_filt(i,j,:)+exp(cmplx(0.,1.)*2*PI*((REAL((i-1)*n)/REAL(nlon-1))+(REAL((j-1)*k)/REAL(nlat-1))))&
!SK201904!!                              *fx*fy*SINC(PI*real(k)/real(a))*SINC(PI*REAL(n)/REAL(a))*BESSEL_JN(1,2*PI*u)/u
!SK201904!!        end do
!SK201904!!        do n=1,a
!SK201904!!          u=sqrt((fx*n)**2+(fy*k)**2)
!SK201904!!          fourier_filt(i,j,:)=fourier_filt(i,j,:)+exp(cmplx(0.,1.)*2*PI*((REAL((i-1)*n)/REAL(nlon-1))+(REAL((j-1)*k)/REAL(nlat-1))))&
!SK201904!!                              *fx*fy*SINC(PI*real(k)/real(a))*SINC(PI*REAL(n)/REAL(a))*BESSEL_JN(1,2*PI*u)/u
!SK201904!!        end do
!SK201904!!        n=0
!SK201904!!        if(k==0) then
!SK201904!!          u=PI
!SK201904!!          fourier_filt(i,j,:)=fourier_filt(i,j,:)+fx*fy*PI
!SK201904!!        else
!SK201904!!          u=sqrt((fx*n)**2+(fy*k)**2)
!SK201904!!          u= BESSEL_JN(1,2*PI*u)/u
!SK201904!!          fourier_filt(i,j,:)=fourier_filt(i,j,:)+exp(cmplx(0.,1.)*2*PI*(REAL((j-1)*k)/REAL(nlat-1)))*fx*fy*SINC(PI*real(k)/real(a))*u
!SK201904!!        end if
!SK201904!!      end do
!SK201904!!      fourier_filt(i,j,:)=fourier_filt(i,j,:)*fourier_vec(i,j,:)
!SK201904!!    end do
!SK201904!!  end do 
!SK201904!!  do i=1,nbv
!SK201904!!    CALL CFFT2B (LDIM, L, M, fourier_filt(:,:,i), WSAVE, LENSAV, WORK, LENWRK, IER)
!SK201904!!    vec_filt(:,:,i)=REAL(fourier_filt(:,:,i),r_size)
!SK201904!!  end do
!SK201904!  deallocate (WORK)
!SK201904!  deallocate (WSAVE)
!SK201904!  RETURN
!SK201904!  END SUBROUTINE LANCZOS_FILTER_vec_naive
END MODULE filter
