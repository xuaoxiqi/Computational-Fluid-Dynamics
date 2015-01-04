
    program main
    implicit none
    real(8), parameter :: Pi=4.0d0*atan(1.0d0)
    real(8), parameter :: Lx=1.0d0
    integer, parameter :: nH=500
    real(8) :: dx
    integer :: i
    real(8) :: X(0:nH)
    real(8) :: uExact(0:nH)

    dx = dble(Lx)/dble(nH)
    do i=0,nH
        X(i) = dble(i)*dx
    enddo

    do i=0,nH
        uExact(i) = 100.0d0*(1.0d0-dexp(25.0d0*(X(i)-1.0d0)))/(1.0d0-dexp(-25.0d0))
    enddo
    open(unit=01,file="exact.dat",status="unknown")
    do i=0,nH
        write(01,*) X(i), uExact(i)
    enddo
    close(01)

    stop
    end program main
