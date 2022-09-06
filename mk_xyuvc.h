write(fname,'("xyuvc",i7.7)') nstep

open(11,file=fname)

do i = 0,Nr
    ! ! inverse
    do j = Ntheta,1,-1
        write(11,'(20e20.10)') xs(i,j),-ys(i,j),us(i,j),-vs(i,j),c(i,j)
    enddo

    do j = 0,Ntheta
        write(11,'(20e20.10)') xs(i,j),ys(i,j),us(i,j),vs(i,j),c(i,j)
    enddo
    write(11,'()')
enddo

close(11)
write(*,*) 'output  : ',fname
