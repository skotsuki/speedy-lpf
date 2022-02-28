MODULE common_lpf
!=======================================================================
!
! [PURPOSE:] Local Particle Filter (LPF)
!            Model Independent Core Module
!
! [REFERENCES:]
!  [1] Ott et al., 2004: A local ensemble Kalman filter for atmospheric
!    data assimilation. Tellus, 56A, 415-428.
!  [2] Hunt et al., 2007: Efficient Data Assimilation for Spatiotemporal
!    Chaos: A Local Ensemble Transform Kalman Filter. Physica D, 230,
!    112-126.
!
! [HISTORY:]
!  04/20/209 Shunji Kotsuki  Created at RIKEN, Kobe
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mtx

  IMPLICIT NONE

  PUBLIC
!=======================================================================
!  LEKF Model Independent Parameters
!=======================================================================
!  INTEGER,PARAMETER :: nbv=16    ! ensemble size
  !INTEGER,SAVE :: nbv=02
  INTEGER,SAVE :: option_aloc = 0 ! SK 20180821 ; 0:normal 1:adaptive localization

CONTAINS
!=======================================================================
!  Main Subroutine of LETKF Core
!   INPUT
!     nobs             : array size, but only first nobsl elements are used
!     nobsl            : total number of observation assimilated at the point
!     hdxb(nobs,nbv)   : obs operator times fcst ens perturbations
!     rdiag(nobs)      : observation error variance
!     rloc(nobs)       : localization weighting function
!     dep(nobs)        : observation departure (yo-Hxb)
!   OUTPUT
!     trans(nbv,nbv) : transformation matrix
!=======================================================================
SUBROUTINE lpf_core(nobs,nobsl,hdxb,rdiag,rloc,dep,parm_infl,gr_common,type_pfmtx,asis,peff,pvec,pmat)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nobs
  INTEGER,INTENT(IN) :: nobsl
  INTEGER,INTENT(IN) :: type_pfmtx ! 0: Roland et al. (2019;MWR), 1: KK's Monte-Carlo Resampling
  REAL(r_size),INTENT(IN)    :: hdxb(1:nobs,1:nbv)
  REAL(r_size),INTENT(IN)    :: rdiag(1:nobs)
  REAL(r_size),INTENT(IN)    :: rloc(1:nobs)
  REAL(r_size),INTENT(IN)    :: dep(1:nobs)
  REAL(r_size),INTENT(IN)    :: gr_common(1:nbv,1:nbv)
  REAL(r_size),INTENT(INOUT) :: parm_infl
  REAL(r_size) :: hdxb_rinv(nobsl,nbv)
  REAL(r_size) :: eivec(nbv,nbv)
  REAL(r_size) :: eival(nbv)
  REAL(r_size) :: pa(nbv,nbv)
  REAL(r_size) :: work1(nbv,nbv)
  REAL(r_size) :: work2(nbv,nobsl)
  REAL(r_size) :: rho
  REAL(r_size) :: parm(4),sigma_o,gain
  REAL(r_size),PARAMETER :: sigma_b = 0.04d0 !error stdev of parm_infl

  REAL(r_size),INTENT(INOUT) :: asis(1:nbv)
  REAL(r_size),INTENT(OUT)   :: peff
  REAL(r_size),INTENT(OUT)   :: pvec(nbv)      ! LPF when updating only ensemble mean
  REAL(r_size),INTENT(OUT)   :: pmat(nbv,nbv)  ! LPF when resampling (mean and ptb)
  REAL(r_size) :: sqpf, qtmp, swgh
  REAL(r_size) :: qpf(1:nbv), acc(1:nbv), pfwgh(1:nbv), rand(1:nbv), grnd(1:nbv), dep2_Ri(1:nbv), pmat_one(nbv,nbv), pmat_suu(nbv,nbv)

  !==> Penny & Miyoshi 2016 (NPG)
  REAL(r_size),PARAMETER :: epsison_lpf  = 1.0d-10 

  !==> KK's Monte-Carlo Resampling
  INTEGER                :: nmonte

  !==> Marginal Particle Filter ; parameters based on personal communication by KK ; M=040
  INTEGER,PARAMETER :: nsu          = 100
  INTEGER,PARAMETER :: a_sample     = 100
  INTEGER,PARAMETER :: b_sample     =  98
  !REAL(8),PARAMETER :: a_factor     =  0.088d0
  !REAL(8),PARAMETER :: b_factor     = -0.088d0
  REAL(8),PARAMETER :: a_factor     =  0.040d0
  REAL(8),PARAMETER :: b_factor     = -0.040d0
  REAL(8)           :: c_factor, d_factor, sss1, sss2, ppp, qqq


  !==> adaptive resampling amplitudes ;;  Potthast et al. (2019; MWR)
  REAL(r_size) :: resample_factor
  REAL(r_size),PARAMETER :: alpha_smooth = 0.05d0 
  REAL(r_size),PARAMETER :: rho0         = 1.0d0
  REAL(r_size),PARAMETER :: rho1         = 1.4d0
  REAL(r_size),PARAMETER :: ccc0         = 0.02d0
  REAL(r_size),PARAMETER :: ccc1         = 0.20d0


  INTEGER :: i,j,k
! LOGICAL :: exist
  IF(nobsl == 0) THEN 
    pvec(:)   = 1.0d0 / dble( nbv )
    peff      = UNDEF
    pmat(:,:) = 0.0d0
    DO j=1,nbv ; pmat(j,j) = 1.0d0 ; END DO ! no resampling 
    RETURN
  ELSE
!-----------------------------------------------------------------------
!  hdxb Rinv
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nobsl
      hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i) * rloc(i)
    END DO
  END DO
!-----------------------------------------------------------------------    
! Likelihood Computation based on Gauss
!-----------------------------------------------------------------------
  ! dep   :: yo     - Hxf(mean)
  ! hdxb  :: Hxf(i) - Hxf(mean)
  ! rdiag :: err*err (i.e., variance)
  ! rloc  :: 0-1

  !CALL calc_pfwgh_norml(nobs,nobsl,nbv,dep,hdxb,rloc,rdiag,pfwgh)
  CALL calc_pfwgh_kkver(nobs,nobsl,nbv,dep,hdxb,rloc,rdiag,asis,pfwgh)

  swgh     = sum( pfwgh(:) )
  pfwgh(:) = pfwgh(:) / swgh
  asis(:)  = pfwgh(:)

  peff     = 1.0d0  / sum( pfwgh(:)**2.0d0 ) ! effective particle size

  acc(:)   = 0.0d0
  acc(1)   = pfwgh(1)
  DO j=2,nbv
    acc(j) = acc(j-1) + pfwgh(j)
  END DO
  pvec(:)   = pfwgh(:)

!-----------------------------------------------------------------------    
! adaptive inflation (Potthast et al. 2018)
!-----------------------------------------------------------------------
  parm = 0.0d0
  DO i=1,nobsl
    parm(1) = parm(1) + dep(i)*dep(i)/rdiag(i) * rloc(i)
  END DO
  DO j=1,nbv
    DO i=1,nobsl
      parm(2) = parm(2) + hdxb_rinv(i,j) * hdxb(i,j)
    END DO
  END DO
  parm(2) = parm(2) / REAL(nbv-1,r_size)
  parm(3) = SUM(rloc(1:nobsl))
!Pottast2018!  parm(4) = (parm(1)-parm(3))/parm(2) ! diagnosed inflation
!Pottast2018!  print '(3f10.4)', parm_infl, parm(4), ( 1.0d0 - alpha_smooth )*parm_infl + alpha_smooth*parm(4)
!Pottast2018!  parm_infl = ( 1.0d0 - alpha_smooth )*parm_infl + alpha_smooth*parm(4)

  parm(4) = (parm(1)-parm(3))/parm(2) - parm_infl
!  sigma_o = 1.0d0/REAL(nobsl,r_size)/MAXVAL(rloc(1:nobsl))
  sigma_o = 2.0d0/parm(3)*((parm_infl*parm(2)+parm(3))/parm(2))**2
  gain = sigma_b**2 / (sigma_o + sigma_b**2)
!CHECK!  print '(4f10.4)', parm_infl, parm(4), gain, parm_infl + gain * parm(4)
  parm_infl = parm_infl + gain * parm(4)

  !==> calc resampling factor
  IF( parm_infl < rho0 )THEN
    resample_factor = ccc0
  ELSE IF( parm_infl > rho1 )THEN
    resample_factor = ccc1
  ELSE
    resample_factor = ccc0 + (ccc1-ccc0)*(parm_infl-rho0)/(rho1-rho0)
  END IF
!CHECK!  print '(a, 4f10.4)', "  resample_factor :: ", rho0, parm_infl, rho1, resample_factor 

!-----------------------------------------------------------------------    
! Resampling with random numbers
!-----------------------------------------------------------------------
  IF( type_pfmtx==0 )THEN ! 0: Roland et al. (2019;MWR), adaptive resampling
    CALL get_resampling_mtx('MR','ON',nbv,acc,pmat)
    !==> inflation by adding noise
    DO j=1,nbv
      grnd(1:nbv) = gr_common(1:nbv,j) * resample_factor
      swgh        = sum( grnd(1:nbv) ) / dble( nbv )
      grnd(1:nbv) = grnd(1:nbv) - swgh                ! so that mean is zero
      pmat(1:nbv,j) =      pmat(1:nbv,j) + grnd(1:nbv)
    END  DO
  ELSE IF( type_pfmtx==1 )THEN ! 1: KK's Monte-Carlo Resampling
    pmat(:,:) = 0.0d0
    !nmonte = nbv    ! 20200528 based on discussion w/ KK 
    nmonte = nbv*5 ! SK's basic test (appendix)
    DO j=1,nmonte ! MONTE-Carlo Resampling M times
      CALL get_resampling_mtx('MR','ON',nbv,acc,pmat_one)
      pmat(:,:) = pmat(:,:) + pmat_one(:,:) / dble(nmonte)
    END DO
  ELSE IF( type_pfmtx==2 )THEN ! 2: Marginal Particle Filter
    !==> get coefficients of MPF
    sss1 = 0.0d0  ;  sss2 = 0.0d0
    do i=1,a_sample
      sss1 = sss1 + a_factor
      sss2 = sss2 + a_factor**2.0d0  ;  end do
    do i=1,b_sample
      sss1 = sss1 + b_factor
      sss2 = sss2 + b_factor**2  ;  end do

    ppp = 1.0d0 - sss1  ;  qqq = 1.0d0 - sss2
    c_factor = 0.5d0 * ( ppp - dsqrt(  ppp**2.0d0 - 2.0d0*(ppp**2.0d0 - qqq) ) )
    d_factor = 0.5d0 * ( ppp + dsqrt(  ppp**2.0d0 - 2.0d0*(ppp**2.0d0 - qqq) ) )

    !CHECK!sss1 = sss1 + c_factor        + d_factor 
    !CHECK!sss2 = sss2 + c_factor**2.0d0 + d_factor**2.0d0
    !CHECK!print *, sss1, sss2, a_factor, b_factor, c_factor, d_factor
    
    pmat(:,:) = 0.0d0
    DO j=1,a_sample ! MONTE-Carlo Resampling M times
      CALL get_resampling_mtx('MR','OF',nbv,acc,pmat_one)
      pmat(:,:) = pmat(:,:) + pmat_one(:,:)*a_factor ; END DO
    DO j=1,b_sample ! MONTE-Carlo Resampling M times
      CALL get_resampling_mtx('MR','OF',nbv,acc,pmat_one)
      pmat(:,:) = pmat(:,:) + pmat_one(:,:)*b_factor ; END DO
      
      CALL get_resampling_mtx('MR','OF',nbv,acc,pmat_one)
      pmat(:,:) = pmat(:,:) + pmat_one(:,:)*c_factor

      !!!CALL get_resampling_mtx('SU','OF',nbv,acc,pmat_one)
      !!!!CALL get_resampling_mtx('MR','OF',nbv,acc,pmat_one)
    pmat_suu(:,:) = 0.0d0
    DO j=1,nsu
      CALL get_resampling_mtx('SU','OF',nbv,acc,pmat_one)
      pmat_suu(:,:) = pmat_suu(:,:) + pmat_one(:,:) / dble(nsu) ; END DO
      
      pmat(:,:) = pmat(:,:) + pmat_suu(:,:)*d_factor

    !!!print '(a)', "========================ave time ============================="
    !!!do i=1,nbv
    !!!  print '(I3,x,40i1,x,40i1,x,40i1,x,5f8.3)', i, int(pmat_one(i,:)+0.00001), int(pmat_suu(i,:)+0.1), int((pmat(i,:)+0.1)), &
    !!!        sum(pmat(:,i)),     sum(pmat(i,:)), sum(pmat_suu(i,:)), pfwgh(i)*dble(nbv), acc(i)
    !!!end do
    !!!print '(7f)', sss1, sss2, a_factor, b_factor, c_factor, d_factor, a_factor*a_sample+b_factor*b_sample+c_factor+d_factor
    !!!stop
  ENDIF
!-----------------------------------------------------------------------
    RETURN
  END IF
END SUBROUTINE lpf_core

!-----------------------------------------------------------------------    
SUBROUTINE get_resampling_mtx(CC,DG,nbv,acc,pmat)
!-----------------------------------------------------------------------
IMPLICIT NONE
  INTEGER     , INTENT(in)  :: nbv
  CHARACTER(2), INTENT(in)  :: CC, DG
  REAL(r_size), INTENT(in)  :: acc(1:nbv)
  REAL(r_size), INTENT(out) :: pmat(nbv,nbv)

  INTEGER      :: i, j, k, init(nbv), inum(nbv)
  REAL(r_size) :: rand(nbv), temp

  DO j=1,nbv ; init(j) = j  ; END DO
  !CALL com_rand(nbv,rand) ! [0-1]

  IF( CC == 'SU' )THEN
    CALL com_rand_seed(nbv,0,rand) ! [0-1]
    temp = rand(1)
    DO i=1,nbv
      rand(i) = ( dble(i) +  temp - 1.0d0 ) / dble(nbv)
    END DO
  ELSE IF( CC == 'MR' )THEN
    CALL com_rand_seed(nbv,0,rand) ! [0-1]
    CALL quick_sort_asnd(rand,init,1,nbv)
  ELSE
    PRINT *, "  ERROR&SROP :: NOT SUCH OPTION in get_resampling_mtx for CC ", CC
    STOP
  END IF

  !==> generate resampling mxm matrix (KK's diagonal-oriented; str)
IF( DG =='ON' )THEN
  inum(:)   = 0
  !(1) :: gen selected particle
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
  !==> generate resampling mxm matrix (non-diagonal-oriented; str, SK's trial)
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

!-----------------------------------------------------------------------    
SUBROUTINE calc_pfwgh_kkver(nobs,nobsl,nbv,dep,hdxb,rloc,rdiag,asis,pfwgh)
!-----------------------------------------------------------------------
IMPLICIT NONE
  INTEGER     , INTENT(IN)    :: nobs,nobsl, nbv
  REAL(r_size), INTENT(IN)    :: dep(1:nobs),hdxb(1:nobs,1:nbv), rloc(1:nobs), rdiag(1:nobs), asis(1:nbv)
  REAL(r_size), INTENT(OUT)   :: pfwgh(1:nbv)
  INTEGER :: i, j, k
  REAL(r_size) :: sqpf, qtmp, qpf(1:nbv), dep2_Ri(1:nbv)

  dep2_Ri(:) = 0.0d0
  DO j=1,nbv
    DO i=1,nobsl
      dep2_Ri(j) = dep2_Ri(j) - 0.5d0*((hdxb(i,j)-dep(i))**2.0d0) * rloc(i) / rdiag(i)
    END DO
  END DO

  sqpf     = 0.0d0
  qpf(:)   = 0.0d0
  DO j=1,nbv
    qtmp = 0.0d0
    DO k=1,nbv
      qtmp = qtmp + dexp( dep2_Ri(k) - dep2_Ri(j) ) * asis(k)
      !print '(2i,5f)', j, k, dep2_Ri(j), dep2_Ri(k), dexp( dep2_Ri(j) ), qtmp, dexp( dep2_Ri(k) - dep2_Ri(j) )
    END DO
    qpf(j) = asis(j) / qtmp
    sqpf   = sqpf + qpf(j)
    !print *, j, qpf(j), dexp( dep2_Ri(j) )
  END DO
  pfwgh(:) = qpf(:) / sqpf
  RETURN
END SUBROUTINE calc_pfwgh_kkver

!-----------------------------------------------------------------------    
SUBROUTINE calc_pfwgh_norml(nobs,nobsl,nbv,dep,hdxb,rloc,rdiag,pfwgh)
!-----------------------------------------------------------------------
IMPLICIT NONE
  INTEGER     , INTENT(IN)  :: nobs,nobsl, nbv
  REAL(r_size), INTENT(IN)  :: dep(1:nobs),hdxb(1:nobs,1:nbv), rloc(1:nobs), rdiag(1:nobs)
  REAL(r_size), INTENT(OUT) :: pfwgh(1:nbv) 
  INTEGER :: i, j, k
  REAL(r_size) :: sqpf, qtmp, qpf(1:nbv), dep2_Ri(1:nbv)

  dep2_Ri(:) = 0.0d0
  DO j=1,nbv
    DO i=1,nobsl
      dep2_Ri(j) = dep2_Ri(j) - 0.5d0*((hdxb(i,j)-dep(i))**2.0d0) * rloc(i) / rdiag(i)
    END DO
  END DO

  sqpf     = 0.0d0
  qpf(:)   = 0.0d0
  DO j=1,nbv
    qpf(j) = dexp( dep2_Ri(j) )
    sqpf   = sqpf + qpf(j)
    !print *, j, qpf(j), dexp( dep2_Ri(j) )
  END DO
  pfwgh(:) = qpf(:) / sqpf
  RETURN
END SUBROUTINE calc_pfwgh_norml
END MODULE common_lpf
