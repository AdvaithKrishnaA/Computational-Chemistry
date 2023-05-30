module q
    real*8, parameter:: pi = acos(-1.d0)
    complex*16, parameter:: iota = (0.d0,1.d0)
    real*8:: alpha, x0, k0, mass, xmin, tmin, dx, dt, dk, pn, pd, tp
    integer:: nx, nt, f
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
end program main

subroutine init()
    use q
    integer:: i,j
    alpha = 20.d0
    x0 = -0.5d0
    k0 = 20.d0
    mass = 14500.d0
    xmin = -2.d0
    tmin = 0.d0
    dx = 0.02d0
    dt = 0.1d0
    nx = 256
    nt = 5000
    pn = 0.d0
    pd = 0.d0
    dk = 2.d0*pi/((nx-1)*dx)
    allocate(x(nx), t(nt), k(nx), V(nx), psi(nx,nt))
    do i = 1, nx
        !space grid
        x(i) = xmin + (i-1)*dx
        !initial wavepacket
        psi(i,1) = ((2*alpha/pi)**0.25)*exp(-alpha*(x(i)-x0)**2)*exp(iota*k0*(x(i)-x0))
        !potential
        if (x(i) < 0.d0 .or. x(i) > 0.5d0) then
            V(i) = 0.d0
        else
            V(i) = 0.1d0
        end if
    end do
    do j = 1, nt
        t(j) = tmin + (j-1)*dt
    end do
    do i = 1, nx
        k(i)= -pi/dx + (i-1)*dk
    end do
    !deallocate(x,t,k,V,psi)
end subroutine init

subroutine tdse()
    use q
    integer:: i, j
    allocate(phi(nx), psi1(nx), psi2(nx), psi3(nx))
    do j = 1, nt-1
        psi1 = psi(:,j)*exp(-iota*V*dt/2.d0)
        call fft(psi1,nx,1)
        psi2 = psi1*exp(-iota*k**2*dt/(2.d0*mass))
        call fft(psi2,nx,-1)
        psi3 = psi2/(real(nx))
        psi(:,j+1) = psi3*exp(-iota*V*dt/2.d0)
    end do
    deallocate(phi, psi1, psi2, psi3)
end subroutine tdse

subroutine forward_dft(sp1, f1, sp2, f2, N)
    integer:: N, j, m
    complex*16:: iota = (0.d0,1.d0)
    real*8, parameter:: pi = acos(-1.d0)
    real*8, dimension(N):: sp1, sp2
    complex*16, dimension(N):: f1, f2
    do m = 1, N
        f2(m) = (0.d0,0.d0)
        do j = 1, N
            f2(m) = f2(m) + f1(j)*exp(-iota*sp2(m)*sp1(j))
        end do
        f2(m) = f2(m)/sqrt(real(N))
    end do
    return
end subroutine forward_dft

subroutine inverse_dft(sp1, f1, sp2, f2, N)
    integer:: N, j, m
    complex*16:: iota = (0.d0,1.d0)
    real*8, parameter:: pi = acos(-1.d0)
    real*8, dimension(N):: sp1, sp2
    complex*16, dimension(N):: f1, f2
    iota = (0.d0,1.d0)
    do j = 1, N
        f2(j) = (0.d0,0.d0)
        do m = 1, N
            f2(j) = f2(j) + f1(m)*exp(iota*sp2(j)*sp1(m))
        end do
        f2(j) = f2(j)/sqrt(real(N))
    end do
    return
end subroutine inverse_dft

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
    tp = pn/pd
    print*, 'transmission probability = ', tp
    f = int(2.d0/dx + 1)
    do i = 1, f
        pn = pn + cdabs(psi(i,5000))*dx
    end do
    tp = pn/pd
    print*, 'reflection probability = ', tp
end subroutine transmission

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