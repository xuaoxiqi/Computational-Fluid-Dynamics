
    program main
    implicit none
    real(8), parameter :: Lx=1.0d0
    real(8) :: dx, dt
    integer :: i
    integer, parameter :: nx=50
    real(8) :: X(0:nx)
    real(8) :: u(0:nx), un(0:nx), uPre(0:nx)
    integer :: itc
    integer, parameter :: itc_max=10000000
    real(8), parameter :: eps=1e-7
    real(8) :: error1, error2, errorU
    real(8), parameter :: constC=0.5d0
    real(8), parameter :: nu=0.02d0

    dx = dble(Lx)/dble(nx)
    dt = 0.01d0
    do i=0,nx
        X(i) = dble(i)*dx
    enddo

    itc = 0
    errorU = 100

    u = 0.0d0
    un = 0.0d0


    do while( (itc.LT.itc_max).AND.(errorU.GT.eps) )

        itc = itc+1

        ! predictor
        do i=1,nx-1
            uPre(i) = u(i)-constC*dt/dx*(u(i+1)-u(i))+nu*dt/dx/dx*(u(i+1)-2.0d0*u(i)+u(i-1))
        enddo
        uPre(0) = 100.0d0
        uPre(nx) = 0.0d0

        ! corrector
        do i=1,nx-1
            u(i) = 0.5d0*( u(i)+uPre(i)-constC*dt/dx*(uPre(i)-uPre(i-1))+nu*dt/dx/dx*(uPre(i+1)-2.0d0*uPre(i)+uPre(i-1)) )
        enddo

        u(0) = 100.0d0
        u(nx) = 0.0d0

        if(MOD(itc,100).EQ.0) then
            error1 = 0.0d0
            error2 = 0.0d0
            do i=1,nx-1
                error1 = error1+dabs(u(i)-un(i))
                error2 = error2+dabs(u(i))
            enddo
            errorU = error1/error2

            un = u
            write(*,*) itc, errorU
        endif

    enddo


    open(unit=01,file="numerical.dat",status="unknown")
    do i=0,nx
        write(01,*) X(i), u(i)
    enddo
    close(01)

    stop
    end program main
