! calculate nc
! nc = c + dt * (adv + vis)


! dc/dr = ((alphainv*c(i+1,j)+(alpha+alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv)
! ddc/drr = ((alphainv*c(i+1,j)+(alpha+alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv)


subroutine cal_nc(Nr,Ntheta,peinv,time,dt,drinvs,ddrinvs,dthetainv,ddthetainv,alpha,alphainv,ur,w,rs,thetas,c,nc)
    implicit none
    integer :: Nr,Ntheta
    double precision :: peinv,time,dt
    double precision :: drinvs(0:Nr-1)
    double precision :: ddrinvs(0:Nr-1)
    double precision :: dthetainv,ddthetainv
    double precision :: alpha,alphainv
    double precision :: ur(0:Nr,0:Ntheta)
    double precision :: w(0:Nr,0:Ntheta)
    double precision :: rs(0:Nr)
    double precision :: thetas(0:Ntheta)
    double precision :: c(0:Nr,0:Ntheta)
    double precision :: nc(0:Nr,0:Ntheta)

    integer :: i,j
    double precision :: adv,vis
    double precision :: advis(0:Nr,0:Ntheta)
    double precision :: rinv,rrinv,palinv
    double precision :: cmax,advismax

    cmax = 1.d-20
    advismax = 1.d-20

    palinv = 1.0d0/(1.0d0+alphainv)

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,rinv,rrinv,adv,vis) &
!$OMP& REDUCTION(max:cmax,advismax)

    do i = 1,Nr-1
        rinv = 1.0d0/rs(i)
        rrinv = 1.0d0/(rs(i)**2)

        j = 0
        adv = -(ur(i,j)*((alphainv*c(i+1,j)+(alpha+alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv))
        vis = (2.0d0*peinv)*( &
            & (2.0d0*alpha*(alphainv*c(i+1,j)-(1.0d0+alphainv)*c(i,j)+c(i-1,j))*ddrinvs(i)*palinv)+ &
            & (2.0d0*rinv*((alphainv*c(i+1,j)+(alpha+alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv)))
        advis(i,j) = adv+vis
        cmax = max(cmax,c(i,j))
        advismax = max(advismax,advis(i,j))

        do j = 1,Ntheta-1
            adv = - &
                & (ur(i,j)* &
                & ((alphainv*c(i+1,j)+(alpha+alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv)+ &
                & w(i,j)*rinv*(c(i,j+1)-c(i,j-1))*0.5d0*dthetainv)

            vis = (2.0d0*peinv)*( &
                & (2.0d0*alpha*(alphainv*c(i+1,j)-(1.0d0+alphainv)*c(i,j)+c(i-1,j))*ddrinvs(i)*palinv)+ &
                & (2.0d0*rinv*((alphainv*c(i+1,j)+(alpha+alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv))+ &
                & (rrinv*(c(i,j+1)-2.0d0*c(i,j)+c(i,j-1))*ddthetainv)+ &
                & ((rrinv/tan(thetas(j)))*(c(i,j+1)-c(i,j-1))*0.5d0*dthetainv))

            advis(i,j) = adv+vis
            cmax = max(cmax,c(i,j))
            advismax = max(advismax,advis(i,j))
        enddo

        j = Ntheta
        adv = -(ur(i,j)*((alphainv*c(i+1,j)+(alpha+alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv))
        vis = (2.0d0*peinv)*( &
            & (2.0d0*alpha*(alphainv*c(i+1,j)-(1.0d0+alphainv)*c(i,j)+c(i-1,j))*ddrinvs(i)*palinv)+ &
            & (2.0d0*rinv*((alphainv*c(i+1,j)+(alpha+alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv)))
        advis(i,j) = adv+vis
        cmax = max(cmax,c(i,j))
        advismax = max(advismax,advis(i,j))
    enddo
!$OMP  END PARALLEL DO

    dt = abs(1.0d0-cmax)*0.9/advismax
    time = time+dt

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)
    do i = 1,Nr-1
        do j = 0,Ntheta
            nc(i,j) = c(i,j)+dt*advis(i,j)
            ! write (*, *) advis(i,j)
        enddo
    enddo
!$OMP  END PARALLEL DO

endsubroutine cal_nc
