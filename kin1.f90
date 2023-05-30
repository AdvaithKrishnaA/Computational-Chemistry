module q
    real*8, parameter:: pi = acos(-1.d0)
    real*8, dimension(:), allocatable:: t, a, b, c, d
    real*8:: dt, k11, k12, k13, k14, k21, k22, k23, k24, k31, k32, k33, k34, k41, k42, k43, k44
    integer:: i, n
end module q

program main
    use q
    implicit none
    dt = 0.0001d0
    n = int(60.d0/dt)
    allocate(t(0:n), a(0:n), b(0:n), c(0:n), d(0:n))
    t(0) = 0.d0
    a(0) = 5.d0
    b(0) = 0.d0
    c(0) = 0.d0
    d(0) = 0.d0
    call rk4()
    deallocate(t, a, b, c, d)
end program main

real*8 function fa(v1, v2)
    real*8:: v1, v2
    fa = -v1 + 3.d0*v2
    return
end function fa

real*8 function fb(v1, v2, v3)
    real*8:: v1, v2, v3
    fb = v1 -3.d0*v2 - 4.2d0*v2 + 7.3d0*v3*v3
    return
end function fb

real*8 function fc(v2, v3)
    real*8:: v2, v3
    fc = 8.4d0*v2 - 14.6d0*v3*v3 - 0.4d0*v3
    return
end function fc

real*8 function fd(v3)
    real*8:: v3
    fd = 0.4d0*v3
    return
end function fd

subroutine rk4()
    use q
    real*8, external:: fa, fb, fc, fd

    open(unit = 11, file = 'data.txt')

    do i = 0, n-1
        t(i+1) = t(i) + dt

        k11 = dt*fa(a(i), b(i))
        k21 = dt*fb(a(i), b(i), c(i))
        k31 = dt*fc(b(i), c(i))
        k41 = dt*fd(c(i))

        k12 = dt*fa(a(i) + k11/2.d0, b(i) + k21/2.d0)
        k22 = dt*fb(a(i) + k11/2.d0, b(i) + k21/2.d0, c(i) + k31/2.d0)
        k32 = dt*fc(b(i) + k21/2.d0, c(i) + k31/2.d0)
        k42 = dt*fd(c(i) + k31/2.d0)

        k13 = dt*fa(a(i) + k12/2.d0, b(i) + k22/2.d0)
        k23 = dt*fb(a(i) + k12/2.d0, b(i) + k22/2.d0, c(i) + k32/2.d0)
        k33 = dt*fc(b(i) + k22/2.d0, c(i) + k32/2.d0)
        k43 = dt*fd(c(i) + k32/2.d0)

        k14 = dt*fa(a(i) + k13, b(i) + k23)
        k24 = dt*fb(a(i) + k13, b(i) + k23, c(i) + k33)
        k34 = dt*fc(b(i) + k23, c(i) + k33)
        k44 = dt*fd(c(i) + k33)

        a(i+1) = a(i) + (k11 + 2.d0*k12 + 2.d0*k13 + k14)/6.d0
        b(i+1) = b(i) + (k21 + 2.d0*k22 + 2.d0*k23 + k24)/6.d0
        c(i+1) = c(i) + (k31 + 2.d0*k32 + 2.d0*k33 + k34)/6.d0
        d(i+1) = d(i) + (k41 + 2.d0*k42 + 2.d0*k43 + k44)/6.d0

        write(11, *) t(i), a(i), b(i), c(i), d(i)
    end do

    close(11)
end subroutine rk4