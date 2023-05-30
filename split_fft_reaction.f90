!Documentation
!Time-dependent wave packet study of non-adiabatic transitions in a A+BC type reaction
!(c) Advaith Krishna A

module q
    real*8, parameter:: pi = acos(-1.d0)
    complex*16, parameter:: iota = (0.d0,1.d0)
    real*8, parameter:: u = 3474.057, xmin = -45, xmax = 45, sigma = 0.3d0, x0 = 9.0d0
    real*8:: dx, dk, dt = 8.d0, p0, b
    integer:: nx = 2048, nt = 1000

    complex*16, dimension(:,:), allocatable:: psi1, psi2
    real*8, dimension(:), allocatable:: x, k, d, V11, V22, V12
    complex*16, dimension(:), allocatable:: phi1, phi2, phi3, phi4
    
    real*8, parameter:: V1 = 4.0167971782296e-2, b1 = 5.5d0, x1 = -4.364721325998e-2, V2 = 4.79833373e-3
    real*8, parameter:: Vasym = 3.61196179e-1, V3= 9.8998917754e-1, b2 = 4.9818195151, x2 = 5.0012635420e-2
    real*8, parameter:: V4 = 1.122019e-2, V5 = 7.9781762366e-1, b3 = 2.3471780470, Vlower = 0.0d0, x3 = -7.6042693477e-1
    real*8, parameter:: b4 = 1.0487590725, x4 = 8.1790045179e-1
end module q

program main
    use q
    implicit none
    allocate(x(nx), k(nx), d(nx), V11(nx), V22(nx), V12(nx), psi1(nx, 0:nt), psi2(nx, 0:nt))
    allocate(phi1(nx), phi2(nx), phi3(nx), phi4(nx))
    call init()
    call evolution()
    call visualise()
    deallocate(x, k, d, V11, V22, V12, psi1, psi2)
    deallocate(phi1, phi2, phi3, phi4)
end program main

subroutine init()
    use q
    integer:: i, flag
    real*8, external:: damp

    dx = (xmax-xmin)/real(nx-1)
    dk = 2.d0*pi/((nx-1)*dx)

    do i = 1, nx
        x(i) = xmin + (i-1)*dx
    end do
    
    call v_init()
    
    flag = int((x0-xmin)/dx + 1)
    p0 = -sqrt(2*u*(0.029 - V11(flag)))
    b = 1./(4.0d0*sigma**2)
    
    do i = 1, nx
        psi1(i, 0) = ((1.d0/(2.d0*pi*sigma**2))**0.25)*exp(-b*(x(i)-x0)**2)*exp(iota*p0*(x(i)-x0))
        psi2(i, 0) = (0.d0, 0.d0)
    end do

    do i = 1, nx/2
        k(i) = dk*(i-1)
    end do

    do i = nx/2 + 1, nx
        k(i) = dk*(i-nx)
    end do
    
    do i = 1, nx
        d(i) = damp(x(i))
    end do

end subroutine init

subroutine v_init()
    use q
    integer:: i
    real*8, external:: f, V1ad, V2ad
    do i = 1, nx
        V11(i) = (1-f(x(i)))*V1ad(x(i)) + f(x(i))*V2ad(x(i))
        V22(i) = f(x(i))*V1ad(x(i)) + (1-f(x(i)))*V2ad(x(i))
        V12(i) = -sqrt(f(x(i))*(1-f(x(i))))*(V2ad(x(i))-V1ad(x(i)))
    end do
end subroutine v_init

subroutine evolution()
    use q
    integer:: j

    do j = 0, nt-1
        phi1=psi1(:,j)
        phi2=psi2(:,j)

        phi1 = phi1*exp(-iota*V11*dt/4.d0)
        phi2 = phi2*exp(-iota*V22*dt/4.d0)

        phi3 = phi1*cos(V12*dt/2.d0) - phi2*iota*sin(V12*dt/2.d0)
        phi4 = -phi1*iota*sin(V12*dt/2.d0) + phi2*cos(V12*dt/2.d0)

        phi1 = phi3*exp(-iota*V11*dt/4.d0)
        phi2 = phi4*exp(-iota*V22*dt/4.d0)

        call fft(phi1,nx,1)
        phi1=phi1*exp(-iota*k*k*dt/(2*u))
        call fft(phi1,nx,-1)
        phi1=phi1/real(nx)

        call fft(phi2,nx,1)
        phi2=phi2*exp(-iota*k*k*dt/(2*u))
        call fft(phi2,nx,-1)
        phi2=phi2/real(nx)

        phi1 = phi1*exp(-iota*V11*dt/4.d0)
        phi2 = phi2*exp(-iota*V22*dt/4.d0)

        phi3 = phi1*cos(V12*dt/2.d0) - phi2*iota*sin(V12*dt/2.d0)
        phi4 = -phi1*iota*sin(V12*dt/2.d0) + phi2*cos(V12*dt/2.d0)

        phi1 = phi3*exp(-iota*V11*dt/4.d0)
        phi2 = phi4*exp(-iota*V22*dt/4.d0)

        psi1(:,j+1) = phi1*d
        psi2(:,j+1) = phi2*d
    end do
    
end subroutine evolution
    

subroutine visualise()
    use q
    integer:: i, j
    character(100):: filename
    real*8, external:: damp

    open(unit = 12, file='data.txt')
    do i = 1, nx
        write(12,*) x(i), V11(i), V22(i), V12(i)
    end do
    close(12)
    
    open(unit = 13, file='psi0.txt')
    do i = 1, nx
        write(13,*) x(i), cdabs(psi1(i,0))**2, cdabs(psi2(i,0))**2
    end do
    close(13)

    open(unit = 14, file='damp.txt')
    do i = 1, nx
        write(14,*) x(i), damp(x(i))
    end do
    close(14)

    do j = 200, 1000, 200
        write(filename,'(i0,a)') j,'.txt'
        open(unit=j, file=filename)
        do i = 1, nx
            write(j,*) x(i), cdabs(psi1(i,j))**2, cdabs(psi2(i,j))**2
        end do
    end do
end subroutine visualise

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

real*8 function f(var)
    use q
    real*8:: var
    f = 0.5d0*(1-tanh(b4*(var - x4)))
    return
end function f

real*8 function V1ad(var)
    use q
    real*8:: var
    V1ad = (V1*exp(b1*(var-x1))/(1 + exp(b1*(var-x1)))**2) + (V2*exp(b1*(var-x1))/(1 + exp(b1*(var-x1))))
    return
end function V1ad

real*8 function V2ad(var)
    use q
    real*8:: var
    V2ad = Vasym - (V3*exp(b2*(var-x2))/(1 + exp(b2*(var-x2)))**2) - (V4*exp(b2*(var-x2))/(1 + exp(b2*(var-x2))))
    V2ad = V2ad - (V5*exp(b3*(var-x3))/(1 + exp(b3*(var-x3)))**2) - Vlower
    return
end function V2ad

real*8 function damp(var)
    use q
    real*8:: var
    if (var .le. -30) then
        damp = sin(pi*(-15 - var)/(30.d0))
    else if (var .ge. 10) then
        damp = sin(pi*(45 - var)/(70.d0))
    else
        damp = 1.d0
    end if
end function damp