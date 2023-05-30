module q
    real*8, parameter:: pi = acos(-1.)
    complex, parameter:: iota = (0.,1.)
    real*8:: xmin, xmax, tmin, tmax, dt, dx
    real*8, dimension(:), allocatable:: x, V
    complex*16, dimension(:,:), allocatable:: psi, k1, k2, k3, k4
    integer:: i, j, n, nx
end module q

program main
    use q
    implicit none
    xmin = 0.
    xmax = 10.
    dx = 0.01
    tmin = 0.
    tmax = 0.03
    dt = 0.0001
    nx = int((xmax-xmin)/dx)
    n = int((tmax-tmin)/dt)
    allocate(x(0:nx), psi(0:nx,0:n), k1(0:nx,0:n), k2(0:nx,0:n), k3(0:nx,0:n), k4(0:nx,0:n), V(0:nx))
    call rk4()
    open(1, file='data.txt')
    write(1,*) 'x t=1 t=150 t=300'
    do i = 0, nx
        write(1,*) x(i), cdabs(psi(i,1))**2, cdabs(psi(i,150))**2, cdabs(psi(i,300))**2
    end do
    close(1)
    deallocate(x,psi,k1,k2,k3,k4,V)
end program main

complex*16 function psi0(v1)
    use q
    real*8:: v1
    !initial wave function
    psi0 = (2.*20./pi)**(1./4.)*exp((-(v1-5.)**2)/(2.*0.01))*exp(20.*iota*v1) !update if needed
end function psi0

complex*16 function f(m,v1,v2)
    use q
    integer:: v1, v2
    complex*16, dimension(0:nx,0:n):: m
    !second derivative using finite difference method
    f = (m(v1+1,v2) + m(v1-1,v2) - 2.*m(v1,v2))/(2.*dx**2) 
end function f

subroutine init()
    use q
    complex*16, external::psi0
    !initial conditions
    do i = 0, nx
        x(i) = xmin + i*dx
        psi(i,0) = psi0(x(i))
        print*, x(i), psi(i,0)
            if (x(i) < 6.d0) then
                V(i) = 0.d0
            else
                V(i) = 1.d0
            end if
    end do
    !boundary conditions
    do j = 0, n
        psi(0,j) = (0.,0.)
        psi(nx,j) = (0.,0.)
    end do
end subroutine init

!modified 4th order Runge-Kutta method for solving TDSE
subroutine rk4()
    use q
    implicit none
    complex*16, external:: f
    call init()
    do j = 0, n-1 !for each time step
        do i = 1, nx-1 !for each spatial step, update k1
            k1(i,j) = iota*f(psi,i,j) - V(i)*psi(i,j)
        end do
        do i = 1, nx-1 !update k2
            k2(i,j) = iota*(f(psi,i,j) + dt/2.*f(k1,i,j)) - V(i)*(psi(i,j)+dt/2.*k1(i,j))
        end do
        do i = 1, nx-1 !update k3
            k3(i,j) = iota*(f(psi,i,j) + dt/2.*f(k2,i,j)) - V(i)*(psi(i,j)+dt/2.*k2(i,j))
        end do
        do i = 1, nx-1 !update k4
            k4(i,j) = iota*(f(psi,i,j) + dt*f(k3,i,j)) - V(i)*(psi(i,j)+dt*k3(i,j))
        end do
        do i = 1, nx-1 !update psi
            psi(i,j+1) = psi(i,j) + dt/6.*(k1(i,j) + 2.*k2(i,j) + 2.*k3(i,j) + k4(i,j))
        end do
    end do
end subroutine rk4