! calculate difference value of advis

subroutine diff(Nr,Ntheta,advis,dif)
    implicit none
    integer :: Nr,Ntheta
    double precision :: advis(0:Nr,0:Ntheta)
    double precision :: dif,dt

    integer :: i,j

    dif = 0.0d0

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j) &
!$OMP& REDUCTION(+:dif)
    do i = 1,Nr-1
        do j = 0,Ntheta
            dif = dif+abs(advis(i,j))
        enddo
    enddo
!$OMP  END PARALLEL DO
    dif = dif/(dble(Nr*Ntheta))
endsubroutine diff
