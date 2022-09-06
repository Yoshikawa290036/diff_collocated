! calculate nc
! nc = c + dt * (adv + vis)

! dc/dr = ((alphainv*c(i+1,j)+(alpha-alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv)
! ddc/drr = ((alphainv*c(i+1,j)+(alpha+alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv)

subroutine cal_advis(Nr,Ntheta,peinv,drinvs,ddrinvs,dthetainv,ddthetainv,alpha,alphainv,ur,w,rs,thetas,advis,c)
    implicit none
    integer :: Nr,Ntheta
    double precision :: peinv
    double precision :: drinvs(0:Nr-1)
    double precision :: ddrinvs(0:Nr-1)
    double precision :: dthetainv,ddthetainv
    double precision :: alpha,alphainv
    double precision :: ur(0:Nr,0:Ntheta)
    double precision :: w(0:Nr,0:Ntheta)
    double precision :: rs(0:Nr)
    double precision :: thetas(0:Ntheta)
    double precision :: c(0:Nr,0:Ntheta)

    integer :: i,j
    double precision :: adv,vis
    double precision :: advis(0:Nr,0:Ntheta)
    double precision :: rinv,rrinv,palinv

    palinv = alpha/(1.0d0+alpha)

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,rinv,rrinv,adv,vis)

    do i = 1,Nr-1
        rinv = 1.0d0/rs(i)
        rrinv = 1.0d0/(rs(i)**2)

        j = 0
        adv = -(ur(i,j)*((alphainv*c(i+1,j)+(alpha-alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv))
        vis = (2.0d0*peinv)*( &
            & (2.0d0*alpha*(alphainv*c(i+1,j)-(1.0d0+alphainv)*c(i,j)+c(i-1,j))*ddrinvs(i)*palinv)+ &
            & (2.0d0*rinv*((alphainv*c(i+1,j)+(alpha-alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv)))
        advis(i,j) = adv+vis

        do j = 1,Ntheta-1
            adv = - &
                & (ur(i,j)* &
                & ((alphainv*c(i+1,j)+(alpha-alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv)+ &
                & w(i,j)*rinv*(c(i,j+1)-c(i,j-1))*0.5d0*dthetainv)

            vis = (2.0d0*peinv)*( &
                & (2.0d0*alpha*(alphainv*c(i+1,j)-(1.0d0+alphainv)*c(i,j)+c(i-1,j))*ddrinvs(i)*palinv)+ &
                & (2.0d0*rinv*((alphainv*c(i+1,j)+(alpha-alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv))+ &
                & (rrinv*(c(i,j+1)-2.0d0*c(i,j)+c(i,j-1))*ddthetainv)+ &
                & ((rrinv/tan(thetas(j)))*(c(i,j+1)-c(i,j-1))*0.5d0*dthetainv))

            advis(i,j) = adv+vis

        enddo

        j = Ntheta
        adv = -(ur(i,j)*((alphainv*c(i+1,j)+(alpha-alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv))
        vis = (2.0d0*peinv)*( &
            & (2.0d0*alpha*(alphainv*c(i+1,j)-(1.0d0+alphainv)*c(i,j)+c(i-1,j))*ddrinvs(i)*palinv)+ &
            & (2.0d0*rinv*((alphainv*c(i+1,j)+(alpha-alphainv)*c(i,j)-alpha*c(i-1,j))*drinvs(i)*palinv)))
        advis(i,j) = adv+vis
    enddo
!$OMP  END PARALLEL DO

endsubroutine cal_advis
