module q
    real*8, parameter:: pi = acos(-1.d0)
    complex*16, parameter:: iota = (0.d0,1.d0)
    integer:: nx, nt
    real*8:: x0, k0, mass, xmin, xmax, tmin, tmax, alpha, dt, dx, dk
    real*8, dimension(:), allocatable:: x, t, k, V
    complex*16, dimension(:,:), allocatable:: psi
    complex*16, dimension(:), allocatable:: phi, psi_der2
end module q

program main
    use q
    implicit none
    call init()
    call tdse()
    call export_setup()
    call shut()
end program main

subroutine init()
    use q
    complex*16, external:: psi0
    integer:: i, j
    !initial parameters
    x0 = -0.5d0
    k0 = 20.d0
    mass = 14500.d0
    xmin = -2.d0
    xmax = 4.d0
    tmin = 0.d0
    tmax = 500.d0
    alpha = 20.d0
    dt = 0.1d0
    dx = 0.02d0
    nx = int((xmax-xmin)/dx)
    nt = int((tmax-tmin)/dt)
    dk = 2.d0*pi/(nx*dx)
    allocate(x(0:nx), t(0:nt), k(0:nx), V(0:nx), psi(0:nx,0:nt))
    do i = 0, nx
        !space grid
        x(i) = xmin + (i)*dx
        !initial wavepacket
        psi(i,0) = psi0(x(i))
        if (x(i) < 0.d0) then
            V(i) = 0.d0
        else
            V(i) = 1.d0
        end if
    end do
    !time
    do j = 0, nt
        t(j) = tmin + (j)*dt
    end do
    !momentum grid
    do i = 0, nx
        k(i) = dk*(i-nx/2)
    end do
end subroutine init

subroutine tdse()
    use q
    integer:: i, j, f
    !euler method for first time step
    allocate(phi(0:nx), psi_der2(0:nx))
    ! call forward_dft(x, psi(0:nx,0), k, phi, nx)
    ! do f = 0, nx
    !     phi(f) = -phi(f)*(k(f)**2)
    ! end do
    ! call inverse_dft(k, phi, x, psi_der2, nx)
    ! do i = 0, nx
    !     psi(i,1) = psi(i,0) + 2.d0*iota*dt*(psi_der2(i)/(2.d0*mass) - V(i)*psi(i,0))
    ! end do
    !time evolution
    do j = 0, nt-1
        call forward_dft(x, psi(0:nx,j), k, phi, nx)
        do f = 0, nx
            phi(f) = -phi(f)*(k(f)**2)
        end do
        call inverse_dft(k, phi, x, psi_der2, nx)
        do i = 0, nx
            psi(i,j+1) = exp(-iota*dt*V(i)/2.d0)*exp(iota*dt*(psi_der2(i)/(2.d0*mass)))*exp(-iota*dt*V(i)/2.d0)*psi(i,j)
        end do   
    end do
    deallocate(phi, psi_der2)
end subroutine tdse

subroutine export_setup()
    use q
    integer:: i,j
    character(100):: filename
    do j = 0, nt, 100
        write(filename,'(i0,a)') j,'.txt'
        open(unit=j, file=filename)
        do i = 0, nx
            write(j,*) x(i), cdabs(psi(i,j))**2
        end do
    end do
end subroutine export_setup

complex*16 function psi0(v1)
    use q
    real*8:: v1
    !initial wave packet
    psi0 = ((2*alpha/pi)**0.25)*exp(-alpha*(v1-x0)**2)*exp(iota*k0*(v1-x0))
end function psi0

subroutine forward_dft(sp1, f1, sp2, f2, N)
    integer:: N, j, m
    complex*16:: iota = (0.d0,1.d0)
    real*8, dimension(0:N):: sp1, sp2
    complex*16, dimension(0:N):: f1, f2
    do m = 0, N
        f2(m) = (0.d0,0.d0)
        do j = 0, N
            f2(m) = f2(m) + f1(j)*exp(-iota*sp2(m)*sp1(j))
        end do
        f2(m) = f2(m)/sqrt(real(N))
    end do
    return
end subroutine forward_dft

subroutine inverse_dft(sp1, f1, sp2, f2, N)
    integer:: N, j, m
    complex*16:: iota = (0.d0,1.d0)
    real*8, dimension(0:N):: sp1, sp2
    complex*16, dimension(0:N):: f1, f2
    iota = (0.d0,1.d0)
    do j = 0, N
        f2(j) = (0.d0,0.d0)
        do m = 0, N
            f2(j) = f2(j) + f1(m)*exp(iota*sp2(j)*sp1(m))
        end do
        f2(j) = f2(j)/sqrt(real(N))
    end do
    return
end subroutine inverse_dft

subroutine shut()
    use q
    deallocate(x, t, k, V, psi)
end subroutine shut