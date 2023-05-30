module q
    real*8, parameter:: pi = acos(-1.)
    real*8:: a, b, h, val,k1,k2,k3,k4
    real*8, dimension(:), allocatable:: x,y
    integer:: i, j, n
end module q

program main
    use q
    implicit none
    a = 0.
    b = 0.2
    h = 0.2
    n = 1
    allocate(x(0:n), y(0:n))
    x(0) = a
    x(n) = b
    y(0) = 1.
    !call euler()
    !call rk2()
    call rk4()
    deallocate(x,y)
end program main

real function f(v1, v2)
    use q
    real*8:: v1, v2
    f = -2*v2 + 4 + v1
    return
end function f

subroutine euler()
    use q
    val = 0.
    do i = 0, n-1
        x(i+1) = x(i) + h
        y(i+1) = y(i) + h * f(x(i), y(i))
    end do
    val = y(n)
    print*,'Euler: ',val
end subroutine euler

subroutine rk2()
    use q
    val = 0.
    do i = 0, n-1
        x(i+1) = x(i) + h
        k1 = f(x(i), y(i))
        k2 = f(x(i+1), y(i) + h*k1)
        y(i+1) = y(i) + h/2 * (k1 + k2)
    end do
    val = y(n)
    print*, 'RK2: ', val
end subroutine rk2

subroutine rk4()
    use q
    val = 0.
    do i = 0, n-1
        x(i+1) = x(i) + h
        k1 = h*f(x(i), y(i))
        k2 = h*f(x(i) + h/2, y(i) + k1/2)
        k3 = h*f(x(i) + h/2, y(i) + k2/2)
        k4 = h*f(x(i+1), y(i) + k3)
        y(i+1) = y(i) + (k1 + 2*k2 + 2*k3 + k4)/6.
    end do
    val = y(n)
    print*, 'RK4: ', val
end subroutine rk4
