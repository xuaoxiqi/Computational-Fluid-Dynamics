
    program main
    implicit none
    real(8), parameter :: Pi=4.0d0*atan(1.0d0)
    real(8), parameter :: Lx=1.0d0, Ly=1.0d0
    real(8) :: dx, dy, dt
    integer :: i, j
    integer, parameter :: nx=21, ny=21
    integer, parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1
    real(8) :: X(nx), Y(ny)
    real(8) :: T(nx,ny), Tn(nx,ny)
    real(8), parameter :: alpha=1.0d0
    integer :: itc
    integer, parameter :: itc_max=10000000
    real(8), parameter :: eps=1e-7
    real(8) :: error1, error2, errorT

    dx = dble(Lx)/dble(nx-1)
    dy = dble(Ly)/dble(ny-1)
    dt = 0.1d0*dx**2
    do i=1,nx
        X(i) = (i-1)*dx
    enddo
    do j=1,ny
        Y(j) = (j-1)*dy
    enddo

    itc = 0
    errorT = 100

    T = 0.0d0
    Tn = 0.0d0

    do while( (itc.LT.itc_max).AND.(errorT.GT.eps) )

        itc = itc+1

        do j=2,ny-1
            do i=2,nx-1
                T(i,j) = T(i,j)+dt*alpha*( (T(i+1,j)-2.0d0*T(i,j)+T(i-1,j))/dx/dx+(T(i,j+1)-2.0d0*T(i,j)+T(i,j-1))/dy/dy )
            enddo
        enddo

        ! Left and right side B.C.
        do j=1,ny
            T(1,j) = 0.0d0
            T(nx,j) = 0.0d0
        enddo
        ! Top and bottom side B.C.
        do i=1,nx
            T(i,1) = (4.0d0*T(i,2)-T(i,3))/3.0d0
            T(i,ny) = dsin(Pi*X(i))
        enddo

        if(MOD(itc,100).EQ.0) then
            error1 = 0.0d0
            error2 = 0.0d0
            do j=2,ny-2
                do i=2,nx-1
                    error1 = error1+dabs(T(i,j)-Tn(i,j))
                    error2 = error2+dabs(T(i,j))
                enddo
            enddo
            errorT = error1/error2

            Tn = T
            write(*,*) itc, errorT
        endif

    enddo




    open(unit=01,file="numerical.dat",status="unknown")
    write(01,*) 'TITLE="2D heat conduction"'
    write(01,*) 'VARIABLES="X" "Y" "T" '
    write(01,101) nx, ny
    do j=1,ny
        do i = 1,nx
            write(01,100) x(i), y(j), T(i,j)
        enddo
    enddo

100 format(2x,10(e12.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')

    close(01)


    open(unit=02,file="yHalf-T.dat",status="unknown")
    do i=1,nx
        write(02,*) X(i), T(i,nyHalf)
    enddo
    close(02)

    open(unit=02,file="xHalf-T.dat",status="unknown")
    do j=1,ny
        write(02,*) T(nxHalf,j), Y(j)
    enddo
    close(02)





    stop
    end program main
