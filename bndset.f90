! set boundary condition of C
! C = 1 at r = 1
! C = 0 at r = R_max
! source  : 杉山先生のメモ (2.4)

! input   : Nr, Ntheta
! return  : c, nc

subroutine bndset(Nr,Ntheta,c)
    implicit none

    integer :: Nr,Ntheta
    double precision :: c(0:Nr,0:Ntheta)

    integer :: j

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(j)
    do j = 0,Ntheta
        ! r = 1
        c(0,j) = 1.0d0

        ! r = R_max
        c(Nr,j) = 0.0d0
    enddo
!$OMP  END PARALLEL DO

endsubroutine bndset
