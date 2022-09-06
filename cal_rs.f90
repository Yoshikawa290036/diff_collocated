! calculate rs, drs

!  C(i-1)          C(i)                  C(i+1)
!     #--------------#----------------------#
!  r(i-1)    ^      r(i)        ^         r(i+1)
!          drs(i-1)           drs(i)

subroutine cal_rs(Nr,drmin,alpha,rs,drinvs,ddrinvs)
    implicit none
    integer :: Nr
    double precision :: drmin,alpha
    double precision :: rs(0:Nr)
    double precision :: drinvs(0:Nr-1)
    double precision :: ddrinvs(0:Nr-1)

    integer :: i
    double precision :: alinv

    rs(0) = 1.0d0
    alinv = 1.0d0/(1.0d0-alpha)

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i)
    do i = 0,Nr-1
        drinvs(i) = 1.0d0/(drmin*(alpha**i))
        ddrinvs(i) = 1.0d0/((drmin*(alpha**i))**2)
        rs(i+1) = 1.0d0+drmin*(1.0d0-alpha**dble(i+1))*alinv
    enddo
!$OMP  END PARALLEL DO

endsubroutine cal_rs
