! calculate Sh
! source  : sugiyama memo (2.5)

! input    : drmininv, dtheta, Nr, Ntheta, thetas, c
! return   : Sh

subroutine cal_Sh(drmininv,dtheta,Nr,Ntheta,thetas,c,Sh)
    implicit none

    double precision :: drmininv,dtheta,Sh
    integer :: Nr,Ntheta
    double precision :: thetas(0:Ntheta)
    double precision :: c(0:Nr,0:Ntheta)

    integer ::  j
    Sh = 0.0d0

!$OMP  PARALLEL DO  REDUCTION(+:Sh) &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(none) &
!$OMP& SHARED(Ntheta,dtheta,drmininv) &
!$OMP& SHARED(c,thetas) &
!$OMP& PRIVATE(j)

    ! do j = 0,Ntheta-1
    !     Sh = Sh-((cos(thetas(j))-cos(thetas(j+1)))*((c(1,j)-c(0,j))*drmininv))
    ! enddo


    do j = 0,Ntheta-1
        Sh = Sh-((cos(thetas(j))-cos(thetas(j+1))))
    enddo
!$OMP  END PARALLEL DO

    ! write(*,*) Sh
    ! stop
endsubroutine cal_Sh
