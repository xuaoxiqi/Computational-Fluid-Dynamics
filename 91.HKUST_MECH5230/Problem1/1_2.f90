
    program main
    implicit none
    real(8), parameter :: Pi=4.0d0*atan(1.0d0)
    real(8), parameter :: Lx=1.0d0
    real(8) :: dx, dt
    integer :: i
    integer, parameter :: nx=15
    integer, parameter :: nxHalf=nx/2
    real(8) :: X(0:nx)
    real(8) :: T(0:nx), Texact(0:nx)
    real(8),parameter :: alpha=0.02d0
    real(8),parameter :: constC=100.d0
    real(8), parameter :: time=20.0d0
    integer :: itc
    integer :: itc_max
    real(8), parameter :: eps=1e-13
    real(8) :: error1, error2, errorT
    real(8), parameter :: cfl=0.06d0

    dx = dble(Lx)/dble(nx)
    dt = 0.1d0
    !dt = cfl*dx*dx/alpha
    write(*,*) "dx=",dx,"   ,dt=",dt
    do i=0,nx
        X(i) = dble(i)*dx
    enddo
    write(*,*) "itc=", time/dt
    itc_max =int(time/dt)

    pause

    itc = 0
    errorT = 100.0d0

    T = constC*dsin(Pi*X/Lx)
    T(0) = 0.0d0
    T(nx) = 0.0d0

    do while( (itc.LE.itc_max).AND.(errorT.GT.eps) )

        do i=1,nx-1
            T(i) = T(i)+dt*alpha*( (T(i+1)-2.0d0*T(i)+T(i-1))/dx/dx )
        enddo

        itc = itc+1

        if(MOD(itc,2).EQ.0) then
            error1 = 0.0d0
            error2 = 0.0d0
            do i=1,nx-1
                Texact(i) = constC*dexp(-alpha*Pi**2*itc*dt/Lx/Lx)*dsin(Pi*X(i)/Lx)
                error1 = error1+dabs(T(i)-Texact(i))**2
                error2 = error2+dabs(Texact(i))**2
            enddo
            errorT = error1/error2

            write(*,*) itc, errorT

            open(unit=02,file='errorT.plt',status="unknown",position='append')
            write(02,*) itc*dt, errorT
            close(02)
        endif

    enddo

    write(*,*) "itc=", itc, "   ,errorT=",errorT


    open(unit=01,file="numerical.dat",status="unknown")
    do i=0,nx
        write(01,*) X(i), T(i)
    enddo
    close(01)

    stop
    end program main
