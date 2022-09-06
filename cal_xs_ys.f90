! calculate xs, ys

subroutine cal_xs_ys(Nr, Ntheta, rs, thetas, xs, ys)
    implicit none
    integer :: Nr, Ntheta
    double precision :: rs(0:Nr)
    double precision :: thetas(0:Ntheta)
    double precision :: xs(0:Nr, 0:Ntheta)
    double precision :: ys(0:Nr, 0:Ntheta)

    integer :: i, j

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(PRIVATE) &
!$OMP& SHARED(Nr,Ntheta,rs,thetas,xs,ys)
    do i = 0, Nr
        do j = 0, Ntheta
            xs(i, j) = rs(i)*cos(thetas(j))
            ys(i, j) = rs(i)*sin(thetas(j))
        end do
    end do
!$OMP  END PARALLEL DO

end subroutine cal_xs_ys
