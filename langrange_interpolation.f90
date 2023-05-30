module q
    real, parameter :: pi=acos(-1.0)
    integer::n,i,j
    real, dimension(:), allocatable::x,y
    real:: l,v,poly
end module q

program prog
    use q
    implicit none
    read(*,*) n
    allocate(x(n), y(n))
    do i=1,n
        read(*,*) x(i), y(i)
    enddo
    write(*,*) 'Enter a value of x'
    read(*,*) v
    call langrange
    write(*,*) 'The value of the polynomial at x=',v,' is ',poly
    deallocate(x,y)
end program prog

subroutine langrange()
    use q
    do i=1,n
        l=1.0
        do j=1,n
            if (i/=j) then
                l=l*(v-x(j))/(x(i)-x(j))
            endif
        enddo
        poly=poly+y(i)*l
    enddo
    return
end subroutine langrange