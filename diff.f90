module q
    real, parameter:: pi=acos(-1.0)
    real:: x, h, fds, bds, cds
end module q

program prog
    use q
    implicit none
    write(*,*) 'Enter the data point and step size'
    read(*,*) x, h
    call der()
    write(*,*) 'FDS:', fds
    write(*,*) 'BDS:', bds
    write(*,*) 'CDS:', cds
end program prog

subroutine der()
    use q
    implicit none
    fds = (sin(x+h)-sin(x))/h
    bds = (sin(x)-sin(x-h))/h
    cds = (sin(x+h)-sin(x-h))/(2*h)
    return
end subroutine der