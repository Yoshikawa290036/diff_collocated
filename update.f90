! update c
! nc : next c
! c <- nc

! input    : Nr, Ntheta, nc
! return   : c

subroutine update(Nr,Ntheta,c,advis,dt)
    implicit none
    integer :: Nr,Ntheta

    double precision :: advis(0:Nr,0:Ntheta)
    double precision :: c(0:Nr,0:Ntheta)
    double precision :: dt

    integer :: i,j

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)
    do i = 1,Nr-1
        do j = 0,Ntheta
            c(i,j) = c(i,j)+dt*advis(i,j)
        enddo
    enddo

!$OMP  END PARALLEL DO
endsubroutine update
