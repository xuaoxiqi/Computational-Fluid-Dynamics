
    program main
    implicit none
    real(8), parameter :: Pi=4.0d0*atan(1.0d0)
    integer, parameter :: Lx=15
    integer, parameter :: nx=20
    real(8) :: dx(nx)
    real(8) :: dt
    integer :: i
    real(8) :: X(0:nx)
    real(8) :: T(0:nx), Tn(0:nx)
    real(8), parameter :: alpha=1.0d0
    integer :: itc
    integer, parameter :: itc_max=10000000
    real(8), parameter :: eps=1e-7
    real(8) :: error1, error2, errorT
    real(8) :: dx_min
    real(8) :: time1, time2
    integer :: itc1, itc2

    do i=0,nx
        X(i) = Lx*(1.0d0-dcos(Pi*i/2.0d0/nx))
    enddo
    dx_min = 100.0d0
    do i=1,nx
        dx(i) = X(i)-X(i-1)
        if(dx(i).LT.dx_min) dx_min = dx(i)
    enddo
    !dt = 3*dx_min**2
    dt = 0.006d0
    write(*,*) dt
    pause

    time1 = 2.0d0
    time2 = 10.0d0
    itc1 = int(time1/dt)
    itc2 = int(time2/dt)

    itc = 0
    errorT = 100

    T = 0.0d0
    Tn = 0.0d0
    T(0) = 0.0d0
    T(nx) = X(nx)**2*dexp(-X(nx))

    do while( (itc.LT.itc_max).AND.(errorT.GT.eps) )

        itc = itc+1

        do i=1,nx-1
            T(i) = T(i)+dt*( alpha*2.0*(dx(i)*(T(i+1)-T(i))-dx(i+1)*(T(i)-T(i-1)))/dx(i)/dx(i+1)/(dx(i)+dx(i+1)) &
                             -(X(i)**2-4.0d0*X(i)+2.0d0)*dexp(-X(i)) )
        enddo

        if(itc.Eq.itc1) then
            open(unit=01,file="t2-non.dat",status="unknown")
            do i=0,nx
                write(01,*) X(i), T(i)
            enddo
            close(01)
        endif


        if(itc.Eq.itc2) then
            open(unit=01,file="t10-non.dat",status="unknown")
            do i=0,nx
                write(01,*) X(i), T(i)
            enddo
            close(01)
        endif

        if(MOD(itc,2).EQ.0) then
            error1 = 0.0d0
            error2 = 0.0d0
            do i=1,nx-1
                error1 = error1+dabs(T(i)-Tn(i))
                error2 = error2+dabs(T(i))
            enddo
            errorT = error1/error2

            Tn = T
            write(*,*) itc, errorT
        endif

    enddo

    open(unit=01,file="numerical.dat",status="unknown")
    do i=0,nx
        write(01,*) X(i), T(i)
    enddo
    close(01)

    open(unit=02,file="grid-non.plt",status="unknown")
    do i=0,nx
         write(02,*) X(i), 2
    enddo
    close(02)

    open(unit=02,file="grid-uni.plt",status="unknown")
    do i=0,nx
         write(02,*) i*Lx/dble(nx), 1
    enddo
    close(02)

    stop
    end program main
