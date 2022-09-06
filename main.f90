program main
    implicit none

    integer :: i,j
    integer :: Nr,Ntheta
    integer :: max_step,nstep
    double precision :: PI,Pe,Sh,alpha,alphainv
    double precision :: lambda,dif,acerror
    double precision :: drmin,drmininv,dtheta,dthetainv,ddthetainv,peinv
    double precision :: time,dt,cfl
    double precision,dimension(:),allocatable :: rs,thetas,drinvs,ddrinvs
    double precision,dimension(:,:),allocatable :: ur,w,c,us,vs,xs,ys,advis
    integer :: imkxyuvc,imkstdout

    character(32) fname
    PI = atan(1.0d0)*4.0d0

! ================================ !
    max_step = 10000
    Nr = 512
    Ntheta = 256
    imkxyuvc = 1000
    imkstdout = 1000
! ================================ !
    drmin = 0.001d0
    alpha = 1.02d0
    acerror = 1.d-10
    read(*,*) Pe
    peinv = 1.0d0/Pe
    alphainv = 1.0d0/alpha
    lambda = 3.0d0/4.0d0
    dtheta = PI/dble(Ntheta)
    dthetainv = 1.0d0/dtheta
    ddthetainv = 1.0d0/(dtheta**2)
    drmininv = 1.0d0/drmin
    nstep = 0
    time = 0.0d0
    dif = 0.0
    dt = 0.0d0
    cfl = 0.9
    include'allocate.h'

    call cal_rs(Nr,drmin,alpha,rs,drinvs,ddrinvs)
    call cal_thetas(Ntheta,dtheta,thetas)

    write(*,'("max step      =           ",1i9)') max_step
    write(*,'("imkuvc        =           ",1i9)') imkxyuvc
    write(*,*)
    write(*,'("Nr            =           ",1i9)') Nr
    write(*,'("Ntheta        =           ",1i9)') Ntheta
    write(*,'("Pe            =           ",20e20.10)') Pe
    write(*,'("Rmax          =           ",20e20.10)') rs(Nr)
    write(*,'("drmin         =           ",20e20.10)') drmin
    write(*,'("alpha         =           ",20e20.10)') alpha
    write(*,'("lambda        =           ",20e20.10)') lambda
    write(*,'("ac error      =           ",20e20.10)') acerror
    write(*,*)
    write(*,*)

    call cal_xs_ys(Nr,Ntheta,rs,thetas,xs,ys)

    call cal_vel(Nr,Ntheta,lambda,rs,thetas,ur,w,us,vs)
    call bndset(Nr,Ntheta,c)
    include'mk_xyuv.h'
    include'mk_xyuvc.h'

    write(fname,'("tdtshdif")')
    open(12,file=fname)
    call cal_Sh(drmininv,dtheta,Nr,Ntheta,thetas,c,Sh)
    write(12,'(20e20.10)') time,dt,Sh,dif

    call cal_dt(pe,drmin,dtheta,dt)

    do nstep = 1,max_step
        call cal_advis(Nr,Ntheta,peinv,drinvs,ddrinvs,dthetainv,ddthetainv,alpha,alphainv,ur,w,rs,thetas,advis,c)

        call update(Nr,Ntheta,c,advis,dt)

        if(mod(nstep,imkstdout)==0) then
            time = dt*dble(nstep)
            call diff(Nr,Ntheta,advis,dif)
            call cal_Sh(drmininv,dtheta,Nr,Ntheta,thetas,c,Sh)
            write(12,'(20e20.10)') time,dt,Sh,dif
            write(*,*) '---------------------------------------'
            write(*,'("nstep        ",1i9.9)') nstep
            write(*,'("time         ",20e20.10)') time
            write(*,'("dt           ",20e20.10)') dt
            write(*,'("Sh           ",20e20.10)') Sh
            write(*,'("dif          ",20e20.10)') dif
            write(*,*)
            if(dif<acerror) then
                write(*,*) 'break time step loop : ',nstep
                exit
            endif
        endif
        if(mod(nstep,imkxyuvc)==0) then
            include'mk_xyuvc.h'
        endif
    enddo
    close(12)

    include'mk_end.h'
endprogram main
