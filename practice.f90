module  q
    real, parameter:: pi=acos(-1.0)
    real, dimension(:), allocatable:: x
    real, dimension(:,:), allocatable:: y
    integer:: i,j,n
    real:: val, poly, dx, p
end module q

program prog
    use q
    implicit none
    write(*,*) 'Enter the number of points'
    read(*,*) n
    allocate(x(n),y(n,n))
    write(*,*) 'Enter a value of x'
    read(*,*) val
    call pol()
    write(*,*) 'The value of the polynomial is', poly
    deallocate(x,y)
end program prog

subroutine grid()
    use q
    dx = (10./real(n-1))
    do i=1,n
        x(i) = dx*i
        y(i,1) = x(i)**3
    end do
end subroutine grid

subroutine interpolation()
    use q
    call grid()
    do j=1,n-1
        do i=1,n-j
            y(i,j+1) = (y(i+1,j)-y(i,j))/(x(i+j)-x(i))
        end do
    end do
end subroutine interpolation

subroutine pol()
    use q
    call interpolation()
    poly = y(1,1)
    do i=2,n
        p = 1.0
        do j=1,i-1
            p = p*(val - x(j))
        end do
        poly = poly + y(1,i)*p
    end do
end subroutine pol
