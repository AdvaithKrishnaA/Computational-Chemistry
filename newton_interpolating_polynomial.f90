module q
    real, parameter :: pi=acos(-1.0)
    integer::n,i,j
    real, dimension(:), allocatable::x
    real, dimension(:,:), allocatable::y
    real:: l,v,pol
end module q

program prog
    use q
    implicit none
    read(*,*) n
    allocate(x(n), y(n,n))
    do i=1,n
        read(*,*) x(i), y(i,1)
    enddo
    write(*,*) 'Enter a value of x'
    read(*,*) v
    call poly
    write(*,*) 'The value of the polynomial at x=',v,' is ',pol
    deallocate(x,y)
end program prog

subroutine coef()
    use q
    do j=1,n-1
        do i=1,n-j
            y(i,j+1)=(y(i+1,j)-y(i,j))/(x(i+j)-x(i))
        enddo
    enddo
    return
end subroutine coef

subroutine poly()
    use q
    call coef()
    pol = y(1,1)
    do i=1,n-1
        p=1.0
        do j=1,i
            p=p*(v-x(j))
        enddo
        pol=pol+y(1,i+1)*p
    enddo
    return
end subroutine poly