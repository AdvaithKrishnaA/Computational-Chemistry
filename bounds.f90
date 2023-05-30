module q
character(1):: JOBZ='V', UPLO='U'
real, parameter:: pi = acos(-1.)
integer:: N=1000, LDA=1000, LWORK, INFO,i,j
real*8, dimension(:,:), allocatable:: A
real*8, dimension(:), allocatable:: W, E
real*8, dimension(:), allocatable:: WORK
real*8:: dx, t
end module q

program fdb
    use q
    implicit none
    dx=5./(N-1)
    t=1./(2.*dx**2)
    LWORK = max(1,3*N-1)
    allocate(A(LDA,N), W(N), E(N), WORK(LWORK))
    call matrix()
    call DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

    open(unit=12, file='data.txt')
    do i=1,N
        E(i) = (2*i*pi)**2/(8*25.)
        write(12,*) i, W(i), E(i)
    enddo
    close(12)

    open(unit=13, file='vector.txt')
    do i=1,N
        write(13,*) i, A(i,1)
    enddo
    close(13)

    deallocate(A, W, E, WORK)
end program fdb

integer function kd(v1,v2)
    use q
    integer:: v1, v2
    if (v1==v2) then
        kd = 1
    else
        kd = 0
    end if
    return
end function kd

subroutine matrix()
    use q
    do i=1,LDA
        do j=1,N
            A(i,j) = (2*t)*kd(i,j) -t*kd(i,(j-1)) -t*kd(i,(j+1))
        enddo
    enddo
    return
end subroutine matrix


























