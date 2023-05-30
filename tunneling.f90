module q
    real, parameter :: pi=acos(-1.0) !pi
    complex, parameter :: i = (0, 1) !iota
    integer::n,ne,j,l
    real, dimension(:), allocatable:: x, V
    complex, dimension(:), allocatable:: psi
    real, dimension(:), allocatable:: E, t
    real:: min, max, dx=0.01, de=0.1, m=1., hbar=1., pavg, val=9.0, pmax, pmin
    complex:: k
end module q

program prog
    use q
    implicit none
    write(*,*) 'Enter the grid boundaries'
    read(*,*) min, max !grid boundaries; here, min = 0 and max = 10
    n = int((max-min)/dx) !number of grid intervals
    ne = int(25/de)
    allocate(x(0:n), V(0:n), psi(0:n), E(0:ne), t(0:ne)) !allocating arrays
    call transmission() !calling subroutines
    call prob()
    deallocate(x, V, psi, E, t) !deallocating arrays
end program prog

subroutine grid(energy)
    use q
    real :: energy !dummy variable for energy
    do j=0,n !initializing the potential
        x(j) = min + j*dx
        if (x(j) < 4.0 .or. x(j) > 5.0) then
            V(j) = 0.0
        else
            V(j) = 9.0 !potential barrier of 9 units
        endif
    enddo
    
    k = csqrt((2,0)*m*(energy-V(1)))/hbar !k = sqrt(2m(E-V))/hbar

    psi(0) = (1,0)
    psi(1) = exp(-i*k*x(1))
    do j=2,n !finite difference method
        psi(j) = (2 - ((2*m/(hbar**2))*(energy-V(j-1))*(dx**2)))*psi(j-1) - psi(j-2)
    enddo
    return
end subroutine grid

subroutine transmission() !calculating transmission probability
    use q
    do l=0,ne
        E(l) = 1 + l*de !energy values
        call grid(E(l))
        pmax = maxval(cabs(psi(600:n))) !max and min values of psi
        pmin = minval(cabs(psi(600:n))) 
        pavg = (pmax**2 + pmin**2)/2 !value of pavg
        t(l) = 2./(1+pavg) !transmission probability
    enddo

    open(unit=199, file='data1.txt')
    do j=0,ne
        write(199,*) E(j), t(j) !writing to file
    enddo
    close(199)
    return
end subroutine transmission

subroutine prob() !calculating probability density for E = 9
    use q
    call grid(val)
    open(unit=200, file='data2.txt')
    do j=0,n
        write(200,*) x(j), cabs(psi(j))**2 !writing to file
    enddo
    close(200)
    return
end subroutine prob

