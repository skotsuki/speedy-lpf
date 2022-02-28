MODULE common_time
  USE common
  USE common_mpi

  IMPLICIT NONE

  PUBLIC

!-----------------------------------------------------------------------
  REAL(r_size), save :: rtimer00, rtimer01, rtimer
  REAL(r_size) :: ltimer00, ltimer01, ltimer
  REAL(r_size) :: ptimer00, ptimer01, ptimer 
  integer, parameter :: slot=20
  REAL(r_size),allocatable,save :: rtimer_all(:,:), rtimerl(:)
  REAL(r_size),allocatable,save :: ltimer_all(:,:), ltimerl(:)

CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE set_timer
  IMPLICIT NONE

  allocate(rtimer_all(slot,nprocs), rtimerl(slot))
  rtimer_all(:,:) = 0.d0
  rtimerl(:) = 0.d0

  allocate(ltimer_all(0:slot,nprocs), ltimerl(0:slot))
  ltimer_all(:,:) = 0.d0
  ltimerl(:) = 0.d0

  RETURN
END SUBROUTINE set_timer
!-----------------------------------------------------------------------
SUBROUTINE flash_timer
  IMPLICIT NONE

  deallocate(rtimer_all, rtimerl)
  deallocate(ltimer_all, ltimerl)

  RETURN
END SUBROUTINE flash_timer
!-----------------------------------------------------------------------
END MODULE common_time
