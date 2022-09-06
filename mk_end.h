

open(13,file='xyuvc_end')

do i = 0,Nr
    ! ! inverse
    do j = Ntheta,1,-1
        write(13,'(20e20.10)') xs(i,j),-ys(i,j),us(i,j),-vs(i,j),c(i,j)
    enddo

    do j = 0,Ntheta
        write(13,'(20e20.10)') xs(i,j),ys(i,j),us(i,j),vs(i,j),c(i,j)
    enddo
    write(13,'()')
enddo

close(13)

write(*,*) 'output  : ',fname
