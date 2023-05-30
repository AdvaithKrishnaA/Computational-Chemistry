module q
    real, parameter :: pi=acos(-1.0)
    integer::num=100,order=4,i,k,s
    real, dimension(100,0:4)::j,n,h
    real, dimension(100)::x
    real:: delx, xmin=0.01, xmax=15
end module q

program prog
    use q
    implicit none
    call printf
end program prog

subroutine CalcDelx()
    use q
    delx = (xmax-xmin)/(num-1)
    return
end subroutine CalcDelx

real function factorial(a)
use q
integer::a
if (a==0) then 
    factorial=1
else 
    factorial=1
    do s=1,a
        factorial=factorial*s
    enddo
endif
return
end function factorial

subroutine bessel()
    use q
    call CalcDelx
    do i=1,num
        x(i)=xmin+(i-1)*delx

        if (x(i)<=0.2) then
            do k=0,order
                j(i,k)=(x(i)**k)/(factorial(k))
            enddo
        else if (x(i)>0.2) then
            j(i,0)=sin(x(i))/x(i)
            j(i,1)=sin(x(i))/(x(i)*x(i))-cos(x(i))/x(i)
            do k=2,order
                j(i,k)=(2*k-1)*j(i,k-1)/x(i)-j(i,k-2)
            enddo
        endif

        n(i,0)=-cos(x(i))/x(i)
        n(i,1)=-cos(x(i))/(x(i)*x(i))-sin(x(i))/x(i)
        do k=2,order
            n(i,k)=(2*k-1)*n(i,k-1)/x(i)-n(i,k-2)
        enddo

        do k=0,order
            h(i,k)=sqrt(j(i,k)**2+n(i,k)**2)
        enddo
    enddo
    return
end subroutine bessel

subroutine printf()
    use q
    call bessel
    open(unit=1,file='data.txt')
    do i=1,num
        write(1,*) x(i), j(i,4), n(i,4), h(i,4)
    enddo
end subroutine printf