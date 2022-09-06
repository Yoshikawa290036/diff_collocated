! calculate thetas

subroutine cal_thetas(Ntheta, dtheta, thetas)
    implicit none

    integer :: Ntheta, j
    double precision :: thetas(0:Ntheta)
    double precision :: dtheta

    !$OMP  PARALLEL DO &
    !$OMP& SCHEDULE(static,1) &
    !$OMP& DEFAULT(SHARED) &
    !$OMP& PRIVATE(j)
    do j = 0, Ntheta
        thetas(j) = dtheta*dble(j)
    end do
    !$OMP  END PARALLEL DO

end subroutine cal_thetas
