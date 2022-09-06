! ! subroutine cal_dt(time,dt)
! subroutine cal_dt(Nr,Ntheta,advis,c,time,dt)

!     implicit none
!     integer :: Nr,Ntheta
!     double precision :: advis(0:Nr,0:Ntheta)
!     double precision :: c(0:Nr,0:Ntheta)
!     double precision :: time,dt,tmp,mmm
!     integer :: i,j

!     mmm = 1.d-2

! !$OMP  PARALLEL DO &
! !$OMP& SCHEDULE(static,1) &
! !$OMP& DEFAULT(SHARED) &
! !$OMP& PRIVATE(i,j,tmp) &
! !$OMP& REDUCTION(min:mmm)
!     do i = 1,Nr-1
!         do j = 0,Ntheta
!             tmp = (1.0d0-c(i,j))/advis(i,j)
!             if(tmp>0) then
!                 mmm = min(mmm,tmp)
!             endif
!         enddo
!     enddo
! !$OMP  END PARALLEL DO

!     dt = 0.1*mmm
!     time = time+dt
! endsubroutine cal_dt


subroutine cal_dt(pe,drmin,dtheta,dt)
    implicit none
    double precision :: pe,drmin,dtheta,dt
    dt = (drmin**2)*pe*0.25d0
    ! dt = min(dt,(dtheta**2)*pe/8.0d0)
endsubroutine cal_dt
