!Documentation
!This program solves the time-dependent Schrodinger equation for a particle in a box with a barrier.
!The transmission and reflection probabilities are calculated.
!The program uses the FFT algorithm to solve the Schrodinger equation.
!(c) Advaith Krishna A

module q
    real*8, parameter:: pi = acos(-1.d0)
    complex*16, parameter:: iota = (0.d0,1.d0)
    real*8:: alpha = 20.d0, x0 = -0.5d0, k0 = 40.d0, mass = 14500.d0 !parameters for the initial wavepacket
    real*8:: xmin = -2.d0, tmin = 0.d0, dx = 0.02d0, dt = 0.1d0, dk, pn, pd, prob !parameters for the grid
    integer:: nx = 256, nt = 5000, f
    real*8, dimension(:), allocatable:: x, t, k, V
    complex*16, dimension(:,:), allocatable:: psi
    complex*16, dimension(:), allocatable:: phi, psi1, psi2, psi3
end module q

program main
    use q
    implicit none
    call init()
    call tdse()
    call transmission()
    call reflection()
end program main

subroutine init()
    use q
    integer:: i
    dk = 2.d0*pi/((nx-1)*dx)
    allocate(x(nx), t(nt), k(nx), V(nx), psi(nx,0:nt))
    do i = 1, nx
        !space grid
        x(i) = xmin + (i-1)*dx
        !initial wavepacket
        psi(i,0) = ((2*alpha/pi)**0.25)*exp(-alpha*(x(i)-x0)**2)*exp(iota*k0*(x(i)-x0))
        !potential
        if (x(i) < 0.d0 .or. x(i) > 0.5d0) then
            V(i) = 0.d0
        else
            V(i) = 0.1d0
        end if
    end do
    do i = 1, nx/2
        k(i)= dk*(i-1)
    end do
    do i = nx/2 + 1, nx
        k(i)= dk*(i-nx)
    end do
end subroutine init

subroutine tdse()
    use q
    integer:: j
    allocate(phi(nx), psi1(nx), psi2(nx), psi3(nx))
    do j = 0, nt-1
        psi1 = psi(:,j)*exp(-iota*V*dt/2.)
        call fft(psi1,nx,1)
        psi2 = psi1*exp(-iota*k**2*dt/(2.*mass))
        call fft(psi2,nx,-1)
        psi3 = psi2/(real(nx))
        psi(:,j+1) = psi3*exp(-iota*V*dt/2.)
    end do
    deallocate(phi, psi1, psi2, psi3)
end subroutine tdse

subroutine transmission()
    use q
    integer:: i
    f = int(2.5d0/dx + 1)
    do i = f, nx
        pn = pn + cdabs(psi(i,5000))*dx
    end do
    do i = 1, nx
        pd = pd + cdabs(psi(i,5000))*dx
    end do
    prob = pn/pd
    print*, 'transmission probability = ', prob
end subroutine transmission

subroutine reflection()
    use q
    integer:: i
    f = int(2.d0/dx + 1)
    do i = 1, f
        pn = pn + cdabs(psi(i,5000))*dx
    end do
    prob = pn/pd
    print*, 'reflection probability = ', prob
end subroutine reflection

subroutine shut()
    use q
    deallocate(x, t, k, V, psi)
end subroutine shut

subroutine fft(x,n,isign)
    implicit real*8(a-h,o-z)
    parameter(npts=10,nptx=2*2**npts)
    complex*16 s,v,w,x(n),cstore
!       dimension of cstore should be at the least twice the grid size
    dimension cstore(nptx)
    complex conjg
    data ntbl/0/
!       roots of unity in the first call
    if(n.gt.ntbl)then
    ntbl=n
    pi=4.0*atan(1.00)
    j=1
    icnt=0
10      s=pi*(0.0,1.0)/float(j)
    do 20 k=0,j-1
    icnt=icnt+1
20      cstore(icnt)=exp(s*float(k))
    j=j+j
    if(j.lt.n)go to 10
    end if
!       permutation of x(j)
    j=1
    do 30 i=1,n
    if(i.le.j)then
    v=x(j)
    x(j)=x(i)
    x(i)=v
    end if
    m=n/2
25      continue
    if(j.gt.m)then
    j=j-m
    m=m/2
    if(m.ge.1)go to 25
    else
    j=j+m
    end if
30       continue
!       multiply x(j) and the roots of unity
    j=1
    icnt=0
40       jj=j+j
     do 50 k=1,j
    icnt=icnt+1
    w=cstore(icnt)
    if(isign.lt.0)w=conjg(w)
    do 50 i=k,n,jj
    v=w*x(i+j)
    x(i+j)=x(i)-v
50      x(i)=x(i)+v
    j=jj
    if(j.lt.n) go to 40

    return
end subroutine fft

!(c) Advaith Krishna A (Indian Institute of Technology Guwahati)