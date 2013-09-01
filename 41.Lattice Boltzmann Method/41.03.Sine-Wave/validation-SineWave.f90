
    program validation
    implicit none
    integer, parameter :: nx=129, ny=129
    real(8), parameter :: Pi=3.14159265359d0
    integer :: i, j
    integer :: nxHalf, nyHalf
    integer :: nxQuarter, nyQuarter
    real(8) :: rho1(nx,ny), rho2(nx,ny)
    real(8) :: tau1, tau2, tau
    real(8) :: rho1Avg, rho2Avg
    real(8) :: delta1, delta2
    real(8) :: r1, r2
    real(8) :: D11, D11Num
    real(8) :: dx
    integer :: nmax
    real(8) :: rho1Exact(nx)
    real(8) :: t

    nxHalf = (nx-1)/2+1
    nyHalf = (ny-1)/2+1
    nxQuarter = (nx-1)/4+1
    nyQuarter = (ny-1)/4+1
    nmax = 2000
    dx = 2.0d0*Pi/(nx-1)
    t = nmax*dx
!-------------------------------------------------------
    open(unit=01,file="sineWave-2000.plt",status="old")
    read(01,*)
    read(01,*)
    read(01,*)
    do j=1,ny
        do i=1,nx
            read(01,*) rho1(i,j), rho2(i,j)
        enddo
    enddo
    close(01)
!--------------------------------------------------------
    tau1 = 0.8d0
    tau2 = 0.6d0
    tau = (tau1+tau2)/2.0d0

    r1 = 1.0d0
    r2 = 1.0d0

    rho1Avg = 1.0d0
    rho2Avg = 1.0d0

    delta1 = 0.001d0
    delta2 = 0.001d0

    D11 = 4.0d0*(tau-0.5d0)*dx/(5.0d0+3.0d0*r1)
    do i=1,nx
        rho1Exact(i) = rho1Avg+delta1*dsin((i-1)*dx)*dexp(-D11*t)
    enddo

!--------------------------------------------------------------
    open(unit=02,file="result.plt",status="unknown")
    write(02,*) 'Title="test"'
    write(02,*) 'VARIABLES="X" "Analytical" "Numerical" '
    write(02,101) nx
    do i=1,nx
        write(02,100) (i-1)*dx, rho1(i,nyHalf), rho1Exact(i)
    enddo

100 format(2x,10(e12.6,'      '))
101 format('ZONE',1x,'I=',1x,i5,2x,'DATAPACKING=POINT')
    close(02)
!--------------------------------------------------------------

    D11Num = dlog(delta1/(rho1(nxQuarter,nxHalf)-rho1Avg))/t

    write(*,*) D11, D11Num

    stop
    end program validation
