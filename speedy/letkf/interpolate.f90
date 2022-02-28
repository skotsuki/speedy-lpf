MODULE interpolate
!=======================================================================
!
! [PURPOSE:] Module for INTERPOLATE with SPEEDY
!
! [HISTORY:]
!   07/03/2018 
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
  PUBLIC :: bilin, inv_interp, rand_grid_gen, unif_grid_gen, inv_weight_gen, bllut_0obs, bllut_me , gauss_filter

CONTAINS

  REAL(r_size) function inv_Ln_distance(x,y,n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: x(2), y(2)
  INTEGER, INTENT(IN) ::  n
  REAL(r_size) :: lattemp(nlon)
  lattemp(1) = -87.159d0
  lattemp(2) = -83.479d0
  lattemp(3) = -79.777d0
  lattemp(4) = -76.070d0
  lattemp(5) = -72.362d0
  lattemp(6) = -68.652d0
  lattemp(7) = -64.942d0
  lattemp(8) = -61.232d0
  lattemp(9) = -57.521d0
  lattemp(10)= -53.810d0
  lattemp(11)= -50.099d0
  lattemp(12)= -46.389d0
  lattemp(13)= -42.678d0
  lattemp(14)= -38.967d0
  lattemp(15)= -35.256d0
  lattemp(16)= -31.545d0
  lattemp(17)= -27.833d0
  lattemp(18)= -24.122d0
  lattemp(19)= -20.411d0
  lattemp(20)= -16.700d0
  lattemp(21)= -12.989d0
  lattemp(22)=  -9.278d0
  lattemp(23)=  -5.567d0
  lattemp(24)=  -1.856d0
  lattemp(25)=   1.856d0
  lattemp(26)=   5.567d0
  lattemp(27)=   9.278d0
  lattemp(28)=  12.989d0
  lattemp(29)=  16.700d0
  lattemp(30)=  20.411d0
  lattemp(31)=  24.122d0
  lattemp(32)=  27.833d0
  lattemp(33)=  31.545d0
  lattemp(34)=  35.256d0
  lattemp(35)=  38.967d0
  lattemp(36)=  42.678d0
  lattemp(37)=  46.389d0
  lattemp(38)=  50.099d0
  lattemp(39)=  53.810d0
  lattemp(40)=  57.521d0
  lattemp(41)=  61.232d0
  lattemp(42)=  64.942d0
  lattemp(43)=  68.652d0
  lattemp(44)=  72.362d0
  lattemp(45)=  76.070d0
  lattemp(46)=  79.777d0
  lattemp(47)=  83.479d0
  lattemp(48)=  87.159d0
  inv_Ln_distance=ABS(ACOS(SIND(lattemp(x(2)))*SIND(lattemp(y(2)))+COSD(lattemp(x(2)))*COSD(lattemp(y(2)))*COSD((360.0d0/nlon)*(x(1)-y(1)))))**(-n)
  END FUNCTION inv_Ln_distance

  SUBROUTINE unif_grid_gen(interval,measdim,thindim,gridloc,measloc)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: interval
  INTEGER, INTENT(IN) :: thindim,measdim
  INTEGER, INTENT(OUT) :: gridloc(2,thindim)
  INTEGER, INTENT(OUT) :: measloc(2,measdim)
  INTEGER :: i,j,k,m,l,o
  k=0
  m=0
  measloc=0.0d0
  gridloc=0.0d0
  DO i=1,nlon,interval
    DO j=1,nlat,interval
      k=k+1
      IF(j+interval<=nlat) THEN
        m=m+1
        GRIDLOC(2,m:(m+interval-2))=(/(l, l=(j+1),(j+interval-1))/)
        GRIDLOC(1,m:(m+interval-2))=i
        m=m+interval-2
      ELSEIF((nlat-j)>0) THEN
        m=m+1
        GRIDLOC(2,m:(m+nlat-j-1))=(/(l, l=(j+1),nlat)/)
        GRIDLOC(1,m:(m+nlat-j-1))=i
        m=m+nlat-j-1
      END IF
      MEASLOC(1,k)=i
      MEASLOC(2,k)=j
    END DO
    IF(i+interval<=nlon) THEN
      DO l=1,(interval-1)
        m=m+1
        GRIDLOC(2,m:(m+nlat-1))=(/(o, o=1,nlat)/)
        GRIDLOC(1,m:(m+nlat-1))=i+l
        m=m+nlat-1
      END DO
    ELSEIF((nlon-i)>0) THEN
      DO l=1,(nlon-i)
        m=m+1
        GRIDLOC(2,m:(m+nlat-1))=(/(o, o=1,nlat)/)
        GRIDLOC(1,m:(m+nlat-1))=i+l
        m=m+nlat-1
      END DO
    END IF
  END DO
  RETURN
  END SUBROUTINE unif_grid_gen

  SUBROUTINE rand_grid_gen(interval,measloc,unmeasloc)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: interval
  INTEGER, INTENT(OUT) :: measloc(2,(1+floor(real(nlon-1)/interval))*(1+floor(real(nlat-1)/interval)))
  INTEGER,INTENT(OUT) :: unmeasloc(2,nlon*nlat-(1+floor(real(nlon-1)/interval))*(1+floor(real(nlat-1)/interval)))
  INTEGER ::  pdfgrid(interval,interval), cdf(interval**2), longrid, latgrid, numberlist(nlon,nlat),i,k
  REAL(r_size) :: randnum

  numberlist=1
  k=0
  CALL RANDOM_SEED()
  DO longrid=1,nlon,interval
    DO latgrid=1,nlat,interval
      pdfgrid=0.0d0
      pdfgrid(1:min(interval,nlon-longrid),1:min(interval,nlat-latgrid))=numberlist(longrid:min(longrid+interval-1,nlon),latgrid:min(latgrid+interval-1,nlat))
      pdfgrid=real(pdfgrid,r_size)/sum(real(pdfgrid,r_size))
      cdf(1)=pdfgrid(1,1)
      CALL RANDOM_NUMBER(randnum)
      i=1
      DO WHILE (randnum>cdf(i))
        i=i+1
        cdf(i)=cdf(i-1)+pdfgrid(1+mod(i-1,interval),(i+interval-1-mod(i+interval-1,interval))/interval)
      END DO
      k=k+1
      measloc(:,k)=(/longrid+mod(i-1,interval),latgrid-1+(interval-1-mod(i+interval-1,interval))/interval/)
      numberlist(measloc(1,k),measloc(2,k))=0
    END DO
  END DO

  k=0
  DO longrid=1,nlon,interval
    DO latgrid=1,nlat,interval
      IF(numberlist(longrid,latgrid)==1) THEN
        k=k+1
        unmeasloc(:,k)=(/longrid,latgrid/)
      END IF
    END DO
  END DO

  RETURN
  END SUBROUTINE rand_grid_gen


  SUBROUTINE inv_weight_gen(interval,measdim,thindim,gridloc,measloc,rho,weight_mat)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: interval, rho
  INTEGER, INTENT(IN) :: thindim,measdim
  INTEGER, INTENT(IN) :: gridloc(2,thindim)
  INTEGER, INTENT(IN) ::  measloc(2,measdim)
  REAL(r_size), INTENT(OUT) :: weight_mat(thindim,measdim)
  INTEGER :: meas, unmeas, sum_ones(measdim)
  INTEGER :: meas_coord(2)
  REAL(r_size) ::  scale_vec(thindim)
  sum_ones=1.0d0
  do meas=1,measdim
    meas_coord=(/measloc(1,meas),measloc(2,meas)/)
    do unmeas=1,thindim
      weight_mat(unmeas,meas)=inv_Ln_distance((/gridloc(1,unmeas),gridloc(2,unmeas)/),meas_coord,rho)
    end do
  end do
  scale_vec=matmul(weight_mat,sum_ones)
  do unmeas=1,thindim
    weight_mat(unmeas,:)=weight_mat(unmeas,:)/scale_vec(unmeas)
  end do
  RETURN
  END SUBROUTINE inv_weight_gen

  SUBROUTINE bilin(vec,vec_int,interval)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: vec(nlon,nlat)
  REAL(r_size), INTENT(OUT) :: vec_int(nlon,nlat)
  REAL(r_size) :: rat, bottom_pole_val, scal_disttomeas, scal_disttobot, disttomeas, disttobot, last_val
  INTEGER, intent(in) :: interval
  INTEGER ::   i, j, k, l, extra, extra2, fin_meas_lat, fin_meas_lon
  IF(interval==1) THEN
    vec_int=vec
    RETURN
  END IF
 
!! ----> If j=lat i=lon, then we interpolate all points in the square ((lat,lon),(lat+interval,lon),(lat,lon+inteval),(lat+interval,lon+interval)) until (lat+n*interval)>nlat-interval, same for lon
 
  RAT=(real(1)/real(INTERVAL)**2)
  VEC_INT=0.0d0
  DO j=1,(nlat-interval),interval
    DO l=0,INTERVAL
      DO i=1,(nlon-interval),interval
        DO k=0,INTERVAL
          VEC_INT(i+k,j+l)=RAT*(VEC(i,j)*(INTERVAL-k)*(INTERVAL-l)+VEC(i+interval,j)*k*(INTERVAL-l)+VEC(i,j+interval)*(INTERVAL-k)*l+VEC(i+interval,j+interval)*k*l)
        END DO
      END DO
    END DO
  END DO

!! ----> For the remaining points to the right, we introduce the first point again after the nlon-th point and interpolate until (lat+n*interval)>nlat-interval
  IF(MOD(nlon,interval)/=1) THEN
    i=i
    extra=nlon-i
    RAT=(real(1)/real(INTERVAL*(EXTRA+1)))
    DO j=1,(nlat-interval),interval
      DO l=0,interval
        DO k=0,extra
          VEC_INT(i+k,j+l)=RAT*(VEC(i,j)*(extra+1-k)*(INTERVAL-l)+VEC(1,j)*k*(INTERVAL-l)+VEC(i,j+interval)*(extra+1-k)*l+VEC(1,j+interval)*k*l)
        END DO
      END DO
    END DO
  END IF

!!----> For the remaining points to the below, we introduce we average the two points above until the (lon+n*interval)>nlon-interval

  IF(MOD(nlat,interval)/=1) THEN
    fin_meas_lat=j
    bottom_pole_val=sum(VEC_INT(:,fin_meas_lat))/nlon
    DO l=1,nlat-fin_meas_lat
      disttobot=real(abs(lat(fin_meas_lat+l)),r_size)
      disttomeas=real(abs(lat(fin_meas_lat)-lat(fin_meas_lat+l)),r_size)
      scal_disttobot=disttobot/(disttobot+disttomeas)
      scal_disttomeas=disttomeas/(disttobot+disttomeas)
      DO i=1,(nlon-interval),interval
        DO k=0,interval
          last_val=VEC_INT(i+k,fin_meas_lat)
          VEC_INT(i+k,fin_meas_lat+l)=bottom_pole_val*scal_disttobot+last_val*scal_disttomeas
        END DO
      END DO
    END DO

!!----> For the remaining points in the bottom right corner, we introduce the first point again and  we average the two points above 

    IF(MOD(nlon,interval)/=1) THEN
      fin_meas_lon=i
      fin_meas_lat=j
      DO l=1,nlat-fin_meas_lat
        disttobot=real(abs(lat(fin_meas_lat+l)),r_size)
        disttomeas=real(abs(lat(fin_meas_lat)-lat(fin_meas_lat+l)),r_size)
        scal_disttobot=disttobot/(disttobot+disttomeas)
        scal_disttomeas=disttomeas/(disttobot+disttomeas)
         DO k=1,nlon-fin_meas_lon
           last_val=VEC_INT(fin_meas_lon+k,fin_meas_lat)
           VEC_INT(fin_meas_lon+k,fin_meas_lat+l)=bottom_pole_val*scal_disttobot+last_val*scal_disttomeas
        END DO
      END DO
    END IF
  END IF
  RETURN
  END SUBROUTINE bilin
   
  
  SUBROUTINE inv_interp(vec,gridloc,measloc,weight_mat,interval,measdim,thindim,vec_interp)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN) :: vec(nlon,nlat)
  REAL(r_size), INTENT(OUT) :: vec_interp(nlon,nlat)
  INTEGER, INTENT(IN) :: thindim,measdim
  REAL(r_size), INTENT(IN) :: weight_mat(thindim,measdim)
  INTEGER, INTENT(IN) :: interval
  INTEGER, INTENT(IN) :: gridloc(2,thindim)
  INTEGER, INTENT(IN) :: measloc(2,measdim)
  REAL(r_size) :: meas_vec_val(measdim)
  REAL(r_size) :: vec_interp_work(thindim)
  INTEGER :: meas, unmeas
  do meas=1,measdim
    meas_vec_val(meas)=vec(measloc(1,meas),measloc(2,meas))
    vec_interp(measloc(1,meas),measloc(2,meas))=vec(measloc(1,meas),measloc(2,meas))
  end do
  vec_interp_work=matmul(weight_mat,meas_vec_val)
  do unmeas=1,thindim
    vec_interp(gridloc(1,unmeas),gridloc(2,unmeas))=vec_interp_work(unmeas)
  end do
  RETURN
  END SUBROUTINE inv_interp

!===================================================================================================================
  SUBROUTINE bllut_0obs(nobsr,vec,msk,wix,wiy,fac)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN)    :: nobsr(nlon,nlat), msk(nlon,nlat), wix(nlon,nlat,4), wiy(nlon,nlat,4), fac(nlon,nlat,4)
  REAL(r_size), INTENT(INOUT) :: vec(nlon,nlat)

  INTEGER       ::   i, j, i1, i2, i3, i4, j1, j2, j3, j4
  REAL(r_size)  :: p1, p2, p3, p4, vec_tmp(nlon,nlat)

  !MEMO!
  ! 20190119 include treatment when no obs
  !          include treatment when wix=-999 & wiy=-999

  vec_tmp(:,:) = vec(:,:)
  DO j=1,nlat  ; DO i=1,nlon
    IF( msk(i,j)>0.5 .and. nobsr(i,j)<0.5 )THEN
      i1=int(wix(i,j,1)+0.1) ; j1=int(wiy(i,j,1)+0.1) ; p1=0.d0 ; if(i1>=1 .and. j1>=1) then ; p1=fac(i,j,1) ; else ;i1=1;j1=1; endif
      i2=int(wix(i,j,2)+0.1) ; j2=int(wiy(i,j,2)+0.1) ; p2=0.d0 ; if(i2>=1 .and. j2>=1) then ; p2=fac(i,j,2) ; else ;i2=1;j2=1; endif
      i3=int(wix(i,j,3)+0.1) ; j3=int(wiy(i,j,3)+0.1) ; p3=0.d0 ; if(i3>=1 .and. j3>=1) then ; p3=fac(i,j,3) ; else ;i3=1;j3=1; endif
      i4=int(wix(i,j,4)+0.1) ; j4=int(wiy(i,j,4)+0.1) ; p4=0.d0 ; if(i4>=1 .and. j4>=1) then ; p4=fac(i,j,4) ; else ;i4=1;j4=1; endif
      !print '(a,2i3,a,8i3,5f9.3,a,5i3)', "inp", i, j, " || ",i1, i2, i3, i4, j1, j2, j3, j4, p1, p2, p3, p4, p1+p2+p3+p4, &
      !                                                " || ",int(nobsr(i,j)),int(nobsr(i1,j1)),int(nobsr(i2,j2)),int(nobsr(i3,j3)),int(nobsr(i4,j4))

      if( nobsr(i1,j1)<0.5 ) p1 = 0.0d0 
      if( nobsr(i2,j2)<0.5 ) p2 = 0.0d0 
      if( nobsr(i3,j3)<0.5 ) p3 = 0.0d0 
      if( nobsr(i4,j4)<0.5 ) p4 = 0.0d0
      if( (p1+p2+p3+p4)>0.1d0 ) &
      vec_tmp(i,j) = (p1*vec(i1,j1)+p2*vec(i2,j2)+p3*vec(i3,j3)+p4*vec(i4,j4) ) / (p1+p2+p3+p4)
    END IF
  END DO       ; END DO
  vec(:,:) = vec_tmp(:,:)
 
  END SUBROUTINE bllut_0obs
!===================================================================================================================
  SUBROUTINE bllut_me(vec,msk_me,wix_me,wiy_me,fac_me,out_me)
  IMPLICIT NONE
  REAL(r_size), INTENT(IN)    :: msk_me(nij1), wix_me(nij1,4), wiy_me(nij1,4), fac_me(nij1,4), vec(nlon,nlat)
  REAL(r_size), INTENT(INOUT) :: out_me(nij1)
  INTEGER ::   ij, i1, i2, i3, i4, j1, j2, j3, j4
  REAL(r_size)  :: p1, p2, p3, p4, tmp_me(nij1)

  tmp_me(:) = out_me(:)
  DO ij=1,nij1
    if( msk_me(ij)>0.5 )then
      out_me(ij) = tmp_me(ij)
    else
      i1=int(wix_me(ij,1)+0.1) ; j1=int(wiy_me(ij,1)+0.1) ; p1=0.d0 ; if(i1>=1 .and. j1>=1) then ; p1=fac_me(ij,1) ; else ;i1=1;j1=1; endif
      i2=int(wix_me(ij,2)+0.1) ; j2=int(wiy_me(ij,2)+0.1) ; p2=0.d0 ; if(i2>=1 .and. j2>=1) then ; p2=fac_me(ij,2) ; else ;i2=1;j2=1; endif
      i3=int(wix_me(ij,3)+0.1) ; j3=int(wiy_me(ij,3)+0.1) ; p3=0.d0 ; if(i3>=1 .and. j3>=1) then ; p3=fac_me(ij,3) ; else ;i3=1;j3=1; endif
      i4=int(wix_me(ij,4)+0.1) ; j4=int(wiy_me(ij,4)+0.1) ; p4=0.d0 ; if(i4>=1 .and. j4>=1) then ; p4=fac_me(ij,4) ; else ;i4=1;j4=1; endif
      out_me(ij) = (p1*vec(i1,j1)+p2*vec(i2,j2)+p3*vec(i3,j3)+p4*vec(i4,j4) ) / (p1+p2+p3+p4)
      !!print '(i,5f20.5)', ij, p1, p2, p3, p4, out_me(ij)
    end if
  END DO
  END SUBROUTINE bllut_me

!===================================================================================================================
  SUBROUTINE gauss_filter( wgh_nij2map, varinp, varout )
  IMPLICIT NONE
  REAL(r_size), INTENT(IN)  :: wgh_nij2map(nij1,nlon,nlat)
  REAL(r_size), INTENT(IN)  :: varinp(nlon,nlat)
  REAL(r_size), INTENT(OUT) :: varout(nij1)

  INTEGER       :: i, j, ij
  REAL(r_size)  :: sfnc, swgh

  varout(:) = 0.0d0
  DO ij=1,nij1
    sfnc=0.0d0 ; swgh=0.0d0
    DO j=1,nlat
      IF( maxval(wgh_nij2map(ij,:,j))>0.0d0 )THEN
        DO i=1,nlon
          IF( wgh_nij2map(ij,i,j)>0.0d0 )THEN
            sfnc = sfnc + wgh_nij2map(ij,i,j)*varinp(i,j)
            swgh = swgh + wgh_nij2map(ij,i,j)
          END IF
        END DO
      END IF
    END DO
    varout(ij) = sfnc / swgh
  END DO
 
  END SUBROUTINE gauss_filter
END MODULE interpolate
