
    program main
    implicit none
    real(8), parameter :: Pi=4.0d0*atan(1.0d0)
    real(8), parameter :: Lx=1.0d0, Ly=1.0d0
    integer, parameter :: nx=201, ny=201
    integer, parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1
    real(8) :: T(nx,ny)
    real(8) :: X(nx), Y(ny)
    integer :: i, j
    real(8) :: dx, dy


    dx = dble(Lx)/dble(nx-1)
    dy = dble(Ly)/dble(ny-1)

    do i=1,nx
        X(i) = (i-1)*dx
    enddo
    do j=1,ny
        Y(j) = (j-1)*dy
    enddo

    do j=1,ny
        do i=1,nx
            T(i,j) = dcosh(Pi*Y(j))/dcosh(Pi)*dsin(Pi*X(i))
        enddo
    enddo


    open(unit=01,file="exact.dat",status="unknown")
    write(01,*) 'TITLE="2D heat conduction"'
    write(01,*) 'VARIABLES="X" "Y" "T_exact" '
    write(01,101) nx, ny
    do j=1,ny
        do i = 1,nx
            write(01,100) x(i), y(j), T(i,j)
        enddo
    enddo

100 format(2x,10(e12.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')

    close(01)



    open(unit=02,file="yHalf-Texact.dat",status="unknown")
    do i=1,nx
        write(02,*) X(i), T(i,nyHalf)
    enddo
    close(02)

    open(unit=02,file="xHalf-Texact.dat",status="unknown")
    do j=1,ny
        write(02,*)  T(nxHalf,j), Y(j)
    enddo
    close(02)

    stop
    end program main
