
    program main
    implicit none
    integer, parameter :: Lx=15
    real(8) :: dx, dt
    integer :: i
    integer, parameter :: nx=20
    real(8) :: X(0:nx)
    real(8) :: T(0:nx), Tn(0:nx)
    real(8), parameter :: alpha=1.0d0
    integer :: itc
    integer, parameter :: itc_max=10000000
    real(8), parameter :: eps=1e-7
    real(8) :: error1, error2, errorT
    real(8) :: A(1:nx-1), B(1:nx-1), C(1:nx-1)
    real(8) :: K(1:nx-1)
    real(8) :: time1, time2
    integer :: itc1, itc2

    dx = dble(Lx)/dble(nx)
    dt = 1.0d0

    time1 = 2.0d0
    time2 = 10.0d0
    itc1 = int(time1/dt)
    itc2 = int(time2/dt)
    write(*,*) itc1, itc2
    pause

    do i=0,nx
        X(i) = dble(i)*dx
    enddo


    itc = 0
    errorT = 100

    T = 0.0d0
    Tn = 0.0d0
    T(0) = 0.0d0
    T(nx) = X(nx)**2*dexp(-X(nx))

    !----------------------------------
    A(1) = 0.0d0
    do i=2,nx-1
        A(i) = alpha*dt/2.0d0/dx/dx
    enddo

    do i=1,nx-1
        B(i) = -(1.0d0+alpha*dt/dx/dx)
    enddo

    do i=1,nx-2
        C(i) = alpha*dt/2.0d0/dx/dx
    enddo
    C(nx-1) = 0.0d0
    !----------------------------------


    do while( (itc.LT.itc_max).AND.(errorT.GT.eps) )

        itc = itc+1

        do i=1,nx-1
            K(i) = -T(i)-alpha*dt/2.0d0/dx/dx*(T(i+1)-2.0d0*T(i)+T(i-1))+(X(i)**2-4.0d0*X(i)+2.0d0)*dexp(-X(i))*dt
        enddo
        K(1) = K(1)-T(0)*alpha*dt/2.0d0/dx/dx
        K(nx-1) = K(nx-1)-T(nx)*alpha*dt/2.0d0/dx/dx

        call Thomas(A(1:nx-1),B(1:nx-1),C(1:nx-1),K(1:nx-1),T(1:nx-1),nx-1)

        T(0) = 0.0d0
        T(nx) = X(nx)**2*dexp(-X(nx))

        if(itc.Eq.itc1) then
            open(unit=01,file="t2-Mac.dat",status="unknown")
            do i=0,nx
                write(01,*) X(i), T(i)
            enddo
            close(01)
        endif


        if(itc.Eq.itc2) then
            open(unit=01,file="t10-Mac.dat",status="unknown")
            do i=0,nx
                write(01,*) X(i), T(i)
            enddo
            close(01)
        endif

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


    open(unit=01,file="numerical.dat",status="unknown")
    do i=0,nx
        write(01,*) X(i), T(i)
    enddo
    close(01)

    stop
    end program main



    subroutine Thomas(coeffA,coeffB,coeffC,coeffF,X,n)
    implicit none
    integer :: n, k
    real(8) :: coeffA(n), coeffB(n), coeffC(n), coeffF(n), X(n)
    real(8) :: A(n), B(n), C(n), F(n)
    real(8) :: t

    A = coeffA
    B = coeffB
    C = coeffC
    F = coeffF

    C(1) = C(1)/B(1)
    F(1) = F(1)/B(1)
    do k=2,n
        t = B(k)-C(k-1)*A(k)
        C(k)=C(k)/t
        F(k)=( F(k)-F(k-1)*A(k) )/t
    enddo
    do k=n-1,1,-1
        F(k)=F(k)-C(k)*F(k+1)
    enddo

    X = F

    return
    end subroutine Thomas

