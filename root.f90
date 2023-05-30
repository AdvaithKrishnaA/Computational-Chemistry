module q
    real*8, parameter:: pi = acos(-1.)
    real*8:: x1, x2, x3, h
end module q

program main
    use q
    implicit none
    x1 = -3.0001
    x2 = 0.
    call bisection()
    call newton()
    call secant()
end program main

real function f(var1)
    use q
    real*8:: var1
    f = var1**2 + var1 - 2 !replace function
    return
end function f

real function df(var1)
    use q
    real*8:: var1
    df = 2.*var1 + 1 !replace function
    return
end function df

subroutine bisection()
    use q
    do
        x3 = (x1 + x2)/2.
        if (f(x1)*f(x3) < 0.) then
            x2 = x3
        else
            x1 = x3
        end if
        if (abs(f(x3)) < 0.001) exit
    end do
    print *, 'root using bisection method = ', x3
end subroutine bisection

subroutine newton()
    use q
    do
        x3 = x1 - f(x1)/df(x1)
        if (abs((x3-x1)/x1) < 1.e-5) exit
        x1 = x3
    end do
    print *, 'root using newton-raphson method = ', x3
end subroutine newton

subroutine secant()
    use q
    do
        x3 = x2 - f(x2)*(x2 - x1)/(f(x2) - f(x1))
        x1 = x2
        x2 = x3
        if (abs((x2-x1)/x1) < 1.e-6) exit
    end do
    print *, 'root using secant method = ', x3
end subroutine secant

subroutine false_position()
    use q
    do
        x3 = x2 - f(x2)*(x2 - x1)/(f(x2) - f(x1))
        if (f(x1)*f(x3) < 0.) then
            x2 = x3
        else
            x1 = x3
        end if
        if (abs((x2-x1)/x1) < 1.e-5) exit
    end do
    print *, 'root using false position method = ', x3
end subroutine false_position