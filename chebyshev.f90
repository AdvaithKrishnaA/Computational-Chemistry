module q
    real, parameter :: pi=acos(-1.0)
    integer::n,i,j
    real, dimension(:), allocatable:: x,c,t
    real:: l,v,pol
end module q

program prog
    use q
    implicit none
    read(*,*) n
    allocate(x(0:n), c(0:n), t(0:n))
    write(*,*) 'Enter a value of x'
    read(*,*) v
    call poly
    write(*,*) 'The value of the polynomial at x=',v,' is ',pol
    deallocate(x,c,t)
end program prog

subroutine xn()
    use q
    do i=0,n
        x(i)=cos((2*i+1)*pi/(2*n+2))
    enddo
    return
end subroutine xn

subroutine cn()
    use q
    call xn()
    do i=0,n
        c(i)=0
        do j=0,n
            c(i)=c(i)+exp(x(j))*cos((2*j+1)*i*pi/(2*n+2))
        enddo
        c(i)=c(i)*2./(n+1)
        if (i==0) then
            c(i)=c(i)/2
        endif
    enddo
    return
end subroutine cn

subroutine tn()
    use q
    do i=0,n
        t(i)=cos(i*acos(v))
    enddo
    return
end subroutine tn

subroutine poly()
    use q
    call cn()
    call tn()
    pol=0
    do i=0,n
        pol=pol+c(i)*t(i)
    enddo
    return
end subroutine poly