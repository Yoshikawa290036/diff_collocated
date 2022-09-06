write (fname, '("xyuv")')
open (10, file=fname)

do i = 0, Nr, 4
    ! inverse
    do j = Ntheta - 1, 1, -2
        write (10, '(20e20.10)') xs(i, j), -ys(i, j), us(i, j), -vs(i, j)
    end do

    do j = 0, Ntheta, 2
        write (10, '(20e20.10)') xs(i, j), ys(i, j), us(i, j), vs(i, j)
    end do
    write (10, '()')
end do

close (10)
write (*, *) 'output  : ', fname
