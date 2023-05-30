module q
    real*8, parameter:: pi = acos(-1.d0)
    real*8:: t, m
    real*8, dimension(:), allocatable:: o1, o2, o3
    real*8:: dt, k1, k2, k3, k4
    real*8:: k11, k12, k13, k21, k22, k23, k31, k32, k33, k14, k24, k34
    integer:: i, n
end module q

program main
    use q
    implicit none

    dt = 0.5d0
    k1 = 3.e-12
    k2 = 1.2e-33
    k3 = 5.5e-4
    k4 = 6.9e-16
    n = 2.e8

    allocate(o1(0:n), o2(0:n), o3(0:n))

    o1(0) = 0.d0
    o2(0) = 1.37e16
    o3(0) = 0.d0
    m = (1.37e16)/(0.22d0)
    t = 0.d0
    call rk4()

    deallocate(o1, o2, o3)
end program main

real*8 function fo1(v1, v2, v3)
    use q
    real*8:: v1, v2, v3
    fo1 = 2.d0*k1*v2 - k2*m*v1*v2 + k3*v3 - k4*v1*v3
    return
end function fo1

real*8 function fo2(v1, v2, v3)
    use q
    real*8:: v1, v2, v3
    fo2 = -k1*v2 - k2*m*v1*v2 + k3*v3 + 2.d0*k4*v1*v3
    return
end function fo2

real*8 function fo3(v1, v2, v3)
    use q
    real*8:: v1, v2, v3
    fo3 = k2*m*v1*v2 - k3*v3 - k4*v1*v3
    return
end function fo3

subroutine rk4()
    use q
    real*8, external:: fo1, fo2, fo3
    open(unit=1, file='rk4.txt')

    do i = 0, n-1
        k11 = dt*fo1(o1(i), o2(i), o3(i))
        k21 = dt*fo2(o1(i), o2(i), o3(i))
        k31 = dt*fo3(o1(i), o2(i), o3(i))

        k12 = dt*fo1(o1(i) + k11/2.d0, o2(i) + k21/2.d0, o3(i) + k31/2.d0)
        k22 = dt*fo2(o1(i) + k11/2.d0, o2(i) + k21/2.d0, o3(i) + k31/2.d0)
        k32 = dt*fo3(o1(i) + k11/2.d0, o2(i) + k21/2.d0, o3(i) + k31/2.d0)
        
        k13 = dt*fo1(o1(i) + k12/2.d0, o2(i) + k22/2.d0, o3(i) + k32/2.d0)
        k23 = dt*fo2(o1(i) + k12/2.d0, o2(i) + k22/2.d0, o3(i) + k32/2.d0)
        k33 = dt*fo3(o1(i) + k12/2.d0, o2(i) + k22/2.d0, o3(i) + k32/2.d0)

        k14 = dt*fo1(o1(i) + k13, o2(i) + k23, o3(i) + k33)
        k24 = dt*fo2(o1(i) + k13, o2(i) + k23, o3(i) + k33)
        k34 = dt*fo3(o1(i) + k13, o2(i) + k23, o3(i) + k33)

        o1(i+1) = o1(i) + (k11 + 2.d0*k12 + 2.d0*k13 + k14)/6.d0
        o2(i+1) = o2(i) + (k21 + 2.d0*k22 + 2.d0*k23 + k24)/6.d0
        o3(i+1) = o3(i) + (k31 + 2.d0*k32 + 2.d0*k33 + k34)/6.d0

        if (mod(i, 10000) == 0) then
            write(1,*) t, o1(i), o2(i), o3(i)
            if (o3(i+1) - o3(i) < 1.e-6) exit
        end if

        t = t + dt

    end do
    close(1)

end subroutine rk4

