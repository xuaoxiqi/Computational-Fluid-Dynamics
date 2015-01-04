
    program main
    implicit none
    integer, parameter :: Lx=15
    integer, parameter :: nH=500
    real(8) :: dx
    integer :: i
    real(8) :: X(0:nH)
    real(8) :: Texact(0:nH)

    dx = dble(Lx)/dble(nH)
    do i=0,nH
        X(i) = dble(i)*dx
    enddo

    do i=0,nH
        Texact(i) = X(i)**2*dexp(-X(i))
    enddo
    open(unit=01,file="exact.dat",status="unknown")
    do i=0,nH
        write(01,*) X(i), Texact(i)
    enddo
    close(01)

    stop
    end program main
