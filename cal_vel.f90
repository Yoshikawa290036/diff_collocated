! calculate velocity
! calculate u_r, u_theta

! input   :  Nr, Ntheta, dr, dtheta, rs, psi
! return  :  ur(u_r), w(u_theta), us, vs

subroutine cal_vel(Nr,Ntheta,lambda,rs,thetas,ur,w,us,vs)
    implicit none
    integer :: Nr,Ntheta
    double precision :: lambda
    double precision :: rs(0:Nr)
    double precision :: thetas(0:Ntheta)
    double precision :: ur(0:Nr,0:Ntheta)
    double precision :: w(0:Nr,0:Ntheta)
    double precision :: us(0:Nr,0:Ntheta)
    double precision :: vs(0:Nr,0:Ntheta)
    integer :: i,j
    double precision :: theta,r,rinv,rrrinv
    ! lambda = (3.0d0*kappa + 2.0d0)/(4.0d0*(kappa + 1.0d0))

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(PRIVATE) &
!$OMP& SHARED(Nr,rs,Ntheta,thetas,ur,w,lambda,us,vs)
    do i = 0,Nr
        r = rs(i)
        rinv = 1.0d0/r
        rrrinv = 1.0d0/(r**3)

        do j = 0,Ntheta
            theta = thetas(j)

            ur(i,j) = cos(theta)*(1.0d0-(2.0d0*lambda*rinv)+((2.0d0*lambda-1.0d0)*rrrinv))
            w(i,j) = sin(theta)*(-1.0d0+(lambda*rinv)+((2.0d0*lambda-1.0d0)*0.5d0*rrrinv))
            us(i,j) = ur(i,j)*cos(theta)-r*sin(theta)*w(i,j)
            vs(i,j) = ur(i,j)*sin(theta)+r*cos(theta)*w(i,j)
        enddo
    enddo
!$OMP  END PARALLEL DO

endsubroutine cal_vel
