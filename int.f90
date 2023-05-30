module q
    real*8, parameter:: pi = acos(-1.)
    real*8:: a, b, h, x, val
    integer:: i, j, n
end module q

program main
    use q
    implicit none
    ! print*, "enter the limits a, b and the number of intervals n"
    ! read*, a, b, n
    a = 0.
    b = pi/2.
    n = 6
    h = (b-a)/n
    call rectangle()
    print*, "rectangle rule: ", val
    call trapezoid()
    print*, "trapezoid rule: ", val
    call simpson()
    print*, "simpson rule: ", val
end program main

real function f(var1)
    real*8:: var1
    f = sqrt(sin(var1)) !replace function
    return
end function f

subroutine rectangle()
    use q
    !composite rectangle rule using fortran
    val = 0.
    do i=0,n-1
        x = a + (i)*h
        val = val + f(x)
    end do
    val = val*h
end subroutine rectangle

subroutine trapezoid()
    use q
    !composite trapezoid rule using fortran
    val = 0.
    do i=1,n-1
        x = a + (i)*h
        val = val + f(x) !replace f(x))
    end do
    val = val + 0.5*(f(a) + f(b)) !replace f(x))
    val = val*h
end subroutine trapezoid

subroutine simpson()
    use q
    !composite simpson rule using fortran
    val = 0.
    do i=1,n-1
        x = a + (i)*h
        if(mod(i,2) == 0) then
            val = val + 2*f(x) !replace f(x))
        else
            val = val + 4*f(x) !replace f(x))
        end if
    end do
    val = val + f(a) + f(b) !replace f(x))
    val = val*h/3.
end subroutine simpson