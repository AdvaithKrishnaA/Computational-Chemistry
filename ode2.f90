module q
    real*8, parameter:: pi = acos(-1.)
    real*8:: a, b, h, val, j1, j2, j3, j4, k1, k2, k3, k4, E
    real*8, dimension(:), allocatable:: t,x,p
    ! real*8, dimension(:,:), allocatable:: M
    integer:: i, j, n
end module q

program main
    use q
    implicit none
    a = 0.
    b = 4.*(2*pi)
    h = 0.02*(2*pi)
    n = 200

    allocate(t(0:n), x(0:n), p(0:n))
    t(0) = a
    t(n) = b
    x(0) = 1.
    p(0) = 0.
    call export_setup()
    deallocate(t,x,p)
end program main

real*8 function f(v1, v2, v3)
    use q
    real*8:: v1, v2, v3
    f = -v2
    return
end function f

subroutine euler()
    use q
    real*8, external :: f
    do i = 0, n-1
        !finding x
        t(i+1) = t(i) + h
        x(i+1) = x(i) + h*p(i)
        p(i+1) = p(i) + h*f(t(i), x(i), p(i))
    end do
end subroutine euler

subroutine rk2()
    use q
    real*8, external :: f
    do i = 0, n-1
        !finding x
        t(i+1) = t(i) + h
        j1 = h*f(t(i), x(i), p(i))
        k1 = h*p(i)
        j2 = h*f(t(i+1), x(i) + k1, p(i) + j1)
        k2 = h*(p(i) + j1)
        x(i+1) = x(i) + (k1 + k2)/2.
        p(i+1) = p(i) + (j1 + j2)/2.
    end do
end subroutine rk2

subroutine rk4()
    use q
    real*8, external :: f
    do i = 0, n-1
        !finding x
        t(i+1) = t(i) + h
        j1 = h*f(t(i), x(i), p(i))
        k1 = h*p(i)
        j2 = h*f(t(i) + h/2, x(i) + k1/2, p(i) + j1/2)
        k2 = h*(p(i) + j1/2)
        j3 = h*f(t(i) + h/2, x(i) + k2/2, p(i) + j2/2)
        k3 = h*(p(i) + j2/2)
        j4 = h*f(t(i+1), x(i) + k3, p(i) + j3)
        k4 = h*(p(i) + j3)
        x(i+1) = x(i) + (k1 + 2*k2 + 2*k3 + k4)/6.
        p(i+1) = p(i) + (j1 + 2*j2 + 2*j3 + j4)/6.
    end do
end subroutine rk4

subroutine export_setup()
    use q
    !for euler
    call euler()
    open(1, file='euler.txt')
    do i = 0, n
        write(1,*) t(i)/(2*pi), x(i), p(i), (x(i)**2 + p(i)**2)
    end do
    close(1)
    !for rk2
    call rk2()
    open(2, file='rk2.txt')
    do i = 0, n
        write(2,*) t(i)/(2*pi), x(i), p(i), (x(i)**2 + p(i)**2)
    end do
    close(2)
    !for rk4
    call rk4()
    open(3, file='rk4.txt')
    do i = 0, n
        write(3,*) t(i)/(2*pi), x(i), p(i), (x(i)**2 + p(i)**2)
    end do
    close(3)
end subroutine export_setup

!GNUPlot Script
! set multiplot
! unset key
! set size 0.2, 0.25
! set origin 0.1, 0
! set yrange [-1:1]
! set xrange [-1:1]
! set xlabel 'x'
! set ylabel 'p'
! set title 'rk4'
! plot 'rk4.txt' using 2:3 with lines
! set size 0.2, 0.25
! set origin 0.3, 0
! set yrange [0.97:1.03]
! set xrange [0:1]
! set xlabel 't/T'
! set ylabel '2E'
! plot 'rk4.txt' using 1:4 with lines
! set size 0.2, 0.25
! set origin 0.5, 0
! set xrange [0:1]
! set yrange [-1:1]
! set xlabel 't/T'
! set ylabel 'x'
! plot 'rk4.txt' using 1:2 with lines
! set size 0.2, 0.25
! set origin 0.1, 0.3
! set yrange [-1:1]
! set xrange [-1:1]
! set xlabel 'x'
! set ylabel 'p'
! set title 'rk2'
! plot 'rk2.txt' using 2:3 with lines
! set size 0.2, 0.25
! set origin 0.3, 0.3
! set yrange [0.97:1.03]
! set xrange [0:1]
! set xlabel 't/T'
! set ylabel '2E'
! plot 'rk2.txt' using 1:4 with lines
! set size 0.2, 0.25
! set origin 0.5, 0.3
! set xrange [0:1]
! set yrange [-1:1]
! set xlabel 't/T'
! set ylabel '2E'
! plot 'rk2.txt' using 1:2 with lines
! set size 0.2, 0.25
! set origin 0.1, 0.6
! set yrange [*:*]
! set xrange [*:*]
! set xlabel 'x'
! set ylabel 'p'
! set title 'euler'
! plot 'euler.txt' using 2:3 with lines
! set size 0.2, 0.25
! set origin 0.3, 0.6
! set xlabel 't/T'
! set ylabel '2E'
! plot 'euler.txt' using 1:4 with lines
! set size 0.2, 0.25
! set origin 0.5, 0.6
! set xlabel 't/T'
! set ylabel '2E'
! plot 'euler.txt' using 1:2 with lines
! unset multiplot