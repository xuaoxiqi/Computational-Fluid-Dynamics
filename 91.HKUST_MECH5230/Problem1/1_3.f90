
    program main
    implicit none
    real(8), parameter :: Pi=4.0d0*atan(1.0d0)
    real(8), parameter :: Lx=1.0d0
    real(8) :: dx, dt
    integer :: i
    integer, parameter :: nx=10
    real(8) :: X(0:nx)
    real(8) :: T(0:nx), Tn(0:nx)
    real(8) :: p(0:nx), q(0:nx)
    real(8) :: pn(0:nx-1), qn(1:nx)
    real(8) :: alpha=0.02d0
    real(8) :: constC=100.d0
    real(8), parameter :: time=10.0d0
    integer :: itc
    integer :: itc_max
    real(8), parameter :: eps=1e-7
    real(8) :: error1, error2, errorT

    dx = dble(Lx)/dble(nx)
    dt = 2.0d0*dx*dx/alpha
    write(*,*) "dx=",dx,"   ,dt=",dt
    do i=0,nx
        X(i) = dble(i)*dx
    enddo
    write(*,*) "itc=", time/dt
    itc_max =int(time/dt)

    itc = 0
    errorT = 100.0d0

    T = constC*dsin(Pi*X/Lx)
    Tn = 0.0d0
    T(0) = 0.0d0
    T(nx) = 0.0d0

    p = constC*dsin(Pi*X/Lx)
    q = constC*dsin(Pi*X/Lx)

    do while( (itc.LE.itc_max).AND.(errorT.GT.eps) )

        pn(0) = 0.0d0
        do i=1,nx-1
            pn(i) = ( alpha*dt/dx/dx*pn(i-1)+(1.0d0-alpha*dt/dx/dx)*p(i)+alpha*dt/dx/dx*p(i+1) )/(1.0d0+alpha*dt/dx/dx)
        enddo

        qn(nx) = 0.0d0
        do i=nx-1,1,-1
            qn(i) = ( alpha*dt/dx/dx*qn(i+1)+(1.0d0-alpha*dt/dx/dx)*q(i)+alpha*dt/dx/dx*q(i-1) )/(1.0d0+alpha*dt/dx/dx)
        enddo


        do i=1,nx-1
            T(i) = 0.5d0*(p(i)+q(i))
            p(i) = pn(i)
            q(i) = qn(i)
        enddo

        itc = itc+1

        if(MOD(itc,1).EQ.0) then
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

    write(*,*) "itc=", itc, "   ,errorT=",errorT


    open(unit=01,file="numerical.dat",status="unknown")
    do i=0,nx
        write(01,*) X(i), T(i)
    enddo
    close(01)

    stop
    end program main
