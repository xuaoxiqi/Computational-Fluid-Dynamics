
    program main
    implicit none
    real(8), parameter :: Pi=4.0d0*atan(1.0d0)
    real(8), parameter :: Lx=1.0d0
    integer, parameter :: nH=500
    integer, parameter :: nxHalf=nH/2
    real(8) :: dx
    integer :: i
    real(8) :: alpha=0.02d0
    real(8) :: constC=100.d0
    real(8) :: X(0:nH)
    real(8) :: Texact(0:nH)
    real(8), parameter :: time=20.0d0
    real(8) :: dt
    integer :: itc, itc_max

    dx = dble(Lx)/dble(nH)
    do i=0,nH
        X(i) = dble(i)*dx
    enddo
    dt =0.1d0
    write(*,*) "itc=", time/dt
    itc_max =int(time/dt)


    itc = 0
    do while(itc.LE.itc_max)
        itc = itc+1
        do i=0,nH
            Texact(i) = constC*dexp(-alpha*Pi**2*itc*dt/Lx/Lx)*dsin(Pi*X(i)/Lx)
        enddo
        open(unit=02,file='centerT.plt',status="unknown",position='append')
        write(02,*) itc*dt, Texact(nxHalf)
        close(02)
    enddo

    open(unit=01,file="exact.dat",status="unknown")
    do i=0,nH
        write(01,*) X(i), Texact(i)
    enddo
    close(01)

    stop
    end program main
