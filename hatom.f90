module q
    real*8, parameter:: pi = acos(-1.)
    real*8:: a, b, h, val, j1, j2, j3, j4, k1, k2, k3, k4, emin, emax, de
    real*8, dimension(:), allocatable:: x,y,p
    real*8, dimension(:), allocatable:: E
    real*8, dimension(:,:), allocatable:: M, R
    integer:: i, j, k, n, nE
end module q

program main
    use q
    implicit none
    !given data
    a = 0.0005
    b = 5
    n = 9999
    h = (b-a)/n

    emin = -0.6
    emax = -0.4
    de = 0.01
    nE = int((emax - emin)/de)

    !allocating memory
    allocate(x(0:n), y(0:n), p(0:n), E(0:nE))
    allocate(M(0:n, 0:nE), R(0:n, 0:nE))

    x(0) = a
    x(n) = b
    y(0) = 0.000001
    p(0) = -1000.0

    !plotting radial distribution function for each E
    do j = 0, nE
        E(j) = emin + j*de
        call rk4(E(j))
        do i = 0, n
            M(i,j) = (x(i)*y(i))**2
        end do
    end do

    open(1, file='data.txt')
    2 format (A, 21(' E=',f5.2))
    write(1,2) 'x', E
    do i = 0, n
        write(1,*) x(i), M(i,:)
    end do
    close(1)

    !plotting R(r) vs r for each E
    do j = 0, nE
        E(j) = emin + j*de
        call rk4(E(j))
        do i = 0, n
            R(i,j) = y(i)
        end do
    end do

    open(2, file='data2.txt')
    write(2,2) 'x', E
    do i = 0, n
        write(2,*) x(i), R(i,:)
    end do
    close(2)

    deallocate(x,y,p,E,M,R)
end program main

!defining 2nd order ODE y"=f(x)
real*8 function f(v1, v2, v3, energy)
    use q
    real*8:: v1, v2, v3, energy
    f = -(2*v3/v1 + 2*(energy + 1/v1)*v2)
    return
end function f

!4th order Runge-Kutta method
subroutine rk4(energy)
    use q
    real*8, external :: f
    real*8:: energy
    do i = 0, n-1
        !finding y
        x(i+1) = x(i) + h
        j1 = h*f(x(i), y(i), p(i), energy)
        k1 = h*p(i)
        j2 = h*f(x(i) + h/2, y(i) + k1/2, p(i) + j1/2, energy)
        k2 = h*(p(i) + j1/2)
        j3 = h*f(x(i) + h/2, y(i) + k2/2, p(i) + j2/2, energy)
        k3 = h*(p(i) + j2/2)
        j4 = h*f(x(i+1), y(i) + k3, p(i) + j3, energy)
        k4 = h*(p(i) + j3)
        y(i+1) = y(i) + (k1 + 2*k2 + 2*k3 + k4)/6.
        p(i+1) = p(i) + (j1 + 2*j2 + 2*j3 + j4)/6.
    end do
end subroutine rk4

!GNUPlot Script
! plot for [col=2:22:5] 'data.txt' using 1:col with lines title columnheader
! plot for [col=2:22:5] 'data2.txt' using 1:col with lines title columnheader