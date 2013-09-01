!!!    This program sloves decay of the sine-wave density profile
!!!    Lattice Boltzmann Equation with MRT-LBE model
!!!    Copyright (C) 2013  Ao Xu
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>


    module alldata
        implicit none
        integer, parameter :: nx=129,ny=129
        real(8), parameter :: cs2=1.0d0/3.0d0
        real(8), parameter :: Pi=3.14159265359d0
        integer :: itc, itc_max
        real(8) :: dx, dy, dt
        real(8) :: tau1, tau2, tau
        real(8) :: r1,r2
        real(8) :: eps, error
        real(8) :: u(nx,ny), v(nx,ny)
        real(8) :: up(nx,ny), vp(nx,ny), rho(nx,ny)
        real(8) :: rho1(nx,ny), rho2(nx,ny)
        real(8) :: omega(0:12), f1(0:12,nx,ny), f2(0:12,nx,ny)
        real(8) :: ex(0:12), ey(0:12)
        data ex/0.0d0,1.0d0,0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0, 2.0d0, 0.0d0, -2.0d0, 0.0d0/
        data ey/0.0d0,0.0d0,1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0, 0.0d0, 2.0d0, 0.0d0, -2.0d0/
    end module alldata

!!!     D2Q9 Lattice Vector Properties:
!!!              6   2   5
!!!                \ | /
!!!              3 - 0 - 1
!!!                / | \
!!!              7   4   8

    program main
    use alldata
    implicit none

!!! input initial data
    dx = 2.0d0*Pi/(nx-1)
    dy = 2.0d0*Pi/(ny-1)
    dt = dx
    r1 = 1.0d0
    r2 = 1.0d0
    tau1 = 0.8d0
    tau2 = 0.6d0
    tau = (tau1+tau2)/2.0d0
    itc = 0
    itc_max = INT(2000)
    eps = 1e-5
    error = 100.0d0

!!! set up initial flow field
    call initial()

    do while((error.GT.eps).AND.(itc.LT.itc_max))

!!! collision step
        call collision()

!!! streaming step
        call streaming()

!!! check convergence
        call check()

!!! output preliminary results
        if(MOD(itc,1000).EQ.0) then
            call output()
        endif

    enddo


!!! output data file
    call output()

    write(*,*)
    write(*,*) '************************************************************'
    write(*,*) 'This program sloves Lid Driven Cavity Flow problem using Lattice Boltzmann Method'
    write(*,*) 'Lattice Boltzmann Equation with BGK approximation'
    write(*,*) 'Consider D2Q9 Particle Discrete Velocity model'
    write(*,*) 'nx =',nx,',       ny =',ny
    write(*,*) 'eps =',eps
    write(*,*) 'itc =',itc
    write(*,*) '************************************************************'
    write(*,*)

    stop
    end program main

!!! set up initial flow field
    subroutine initial()
    use alldata
    implicit none
    integer :: i, j
    integer :: alpha
    real(8) :: delta1, delta2
    real(8) :: rho1Avg, rho2Avg
    real(8) :: feq

    omega(0) = 3.0d0/8.0d0
    do alpha=1,4
        omega(alpha) = 1.0d0/12.0d0
    enddo
    do alpha=5,8
        omega(alpha) = 1.0d0/16.0d0
    enddo
    do alpha=9,12
        omega(alpha) = 1.0d0/96.0d0
    enddo


    rho1Avg = 1.0d0
    rho2Avg = 1.0d0
    delta1 = 0.001d0
    delta2 = 0.001d0
    do i=1,nx
        do j=1,ny
            rho1(i,j) = rho1Avg+delta1*dsin(i*dx)
            rho2(i,j) = rho2Avg-delta2*dsin(i*dx)
        enddo
    enddo

    do i=1,nx
        do j=1,ny
            do alpha=0,12
                f1(alpha,i,j) = feq(omega(alpha),alpha,rho1(i,j),r1,u(i,j),v(i,j))
                f2(alpha,i,j) = feq(omega(alpha),alpha,rho2(i,j),r2,u(i,j),v(i,j))
            enddo
        enddo
    enddo

    return
    end subroutine initial

    function feq(omega,k,rho,rk,u,v)
    implicit none
    real(8) :: feq
    integer :: k, alpha
    real(8) :: rk
    real(8) :: u, v, rho
    real(8) :: A, B, C
    real(8) :: us2
    real(8) :: un(0:12)
    real(8) :: omega
    real(8) :: ex(0:12), ey(0:12)
    data ex/0.0d0,1.0d0,0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0, 2.0d0, 0.0d0, -2.0d0, 0.0d0/
    data ey/0.0d0,0.0d0,1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0, 0.0d0, 2.0d0, 0.0d0, -2.0d0/

    alpha = k

    A = 8.0d0/(5.0d0+3.0d0*rk)
    B = (3.0d0*rk-1.0d0)/(5.0d0+3.0d0*rk)
    C = (7.0d0-3.0d0*rk)/(5.0d0+3.0d0*rk)
    us2 = u*u+v*v
    do alpha=0,12
        un(alpha) = u*ex(alpha)+v*ey(alpha)
    enddo

    if(alpha.EQ.0) then
        feq = omega*rho*(A*rk-us2)
    elseif((alpha.GE.1).AND.(alpha.LE.4)) then
        feq = omega*rho*( A+8.0d0*B*un(alpha)+2.0d0*un(alpha)*un(alpha)-us2 )
    elseif((alpha.GE.5).AND.(alpha.LE.8)) then
        feq = omega*rho*( A+2.0d0*A*un(alpha)+2.0d0*un(alpha)*un(alpha)-us2 )
    else
        feq = omega*rho*( A+4.0d0*C*un(alpha)+2.0d0*un(alpha)*un(alpha)-us2 )
    endif

    return
    end function feq

!!! collision step
    subroutine collision()
    use alldata
    implicit none
    integer :: i, j
    integer :: alpha
    real(8) :: g1(0:12,nx,ny), g2(0:12,nx,ny), s1(0:12), s2(0:12)
    real(8) :: geq

    do i=1,nx
        do j=1,ny

            rho1(i,j) = 0.0d0
            do alpha=0,12
                rho1(i,j) = rho1(i,j)+f1(alpha,i,j)
            enddo
            rho2(i,j) = 0.0d0
            do alpha=0,12
                rho2(i,j) = rho2(i,j)+f2(alpha,i,j)
            enddo
            rho(i,j) = rho1(i,j)+rho2(i,j)

            u(i,j) = 0.0d0
            v(i,j) = 0.0d0
            do alpha =0,12
                u(i,j) = u(i,j)+ex(alpha)*(f1(alpha,i,j)+f2(alpha,i,j))
                v(i,j) = v(i,j)+ey(alpha)*(f1(alpha,i,j)+f2(alpha,i,j))
            enddo
            u(i,j) = u(i,j)/rho(i,j)
            v(i,j) = v(i,j)/rho(i,j)

    g1(0,i,j) = f1(0,i,j)+f1(1,i,j)+f1(2,i,j)+f1(3,i,j)+f1(4,i,j)+f1(5,i,j)+f1(6,i,j) &
                +f1(7,i,j)+f1(8,i,j)+f1(9,i,j)+f1(10,i,j)+f1(11,i,j)+f1(12,i,j)
    g1(1,i,j) = f1(1,i,j)-f1(3,i,j)+f1(5,i,j)-f1(6,i,j)-f1(7,i,j)+f1(8,i,j)+2.0d0*f1(9,i,j)-2.0d0*f1(11,i,j)
    g1(2,i,j) = f1(2,i,j)-f1(4,i,j)+f1(5,i,j)+f1(6,i,j)-f1(7,i,j)-f1(8,i,j)+2.0d0*f1(10,i,j)-2.0d0*f1(12,i,j)
    g1(3,i,j) = f1(5,i,j)-f1(6,i,j)+f1(7,i,j)-f1(8,i,j)
    g1(4,i,j) = f1(1,i,j)+f1(2,i,j)+f1(3,i,j)+f1(4,i,j)+2.0d0*( f1(5,i,j)+f1(6,i,j)+f1(7,i,j)+f1(8,i,j) ) &
                +4.0d0*( f1(9,i,j)+f1(10,i,j)+f1(11,i,j)+f1(12,i,j) )
    g1(5,i,j) = f1(1,i,j)-f1(2,i,j)+f1(3,i,j)-f1(4,i,j)+4.0d0*( f1(9,i,j)-f1(10,i,j)+f1(11,i,j)-f1(12,i,j) )
    g1(6,i,j) = f1(1,i,j)-f1(3,i,j)+2.0d0*( f1(5,i,j)-f1(6,i,j)-f1(7,i,j)+f1(8,i,j) )+8.0d0*( f1(9,i,j)-f1(11,i,j))
    g1(7,i,j) = f1(2,i,j)-f1(4,i,j)+2.0d0*( f1(5,i,j)+f1(6,i,j)-f1(7,i,j)-f1(8,i,j) )+8.0d0*( f1(10,i,j)-f1(12,i,j) )
    g1(8,i,j) = f1(1,i,j)-f1(3,i,j)+f1(5,i,j)-f1(6,i,j)-f1(7,i,j)+f1(8,i,j)+8.0d0*( f1(9,i,j)-f1(11,i,j) )
    g1(9,i,j) = f1(2,i,j)-f1(4,i,j)+f1(5,i,j)+f1(6,i,j)-f1(7,i,j)-f1(8,i,j)+8.0d0*( f1(10,i,j)-f1(12,i,j) )
    g1(10,i,j) = f1(1,i,j)+f1(2,i,j)+f1(3,i,j)+f1(4,i,j)+2.0d0*( f1(5,i,j)+f1(6,i,j)+f1(7,i,j)+f1(8,i,j) ) &
                + 16.0d0*( f1(9,i,j)+f1(10,i,j)+f1(11,i,j)+f1(12,i,j) )
    g1(11,i,j) = f1(1,i,j)-f1(2,i,j)+f1(3,i,j)-f1(4,i,j)+16.0d0*( f1(9,i,j)-f1(10,i,j)+f1(11,i,j)-f1(12,i,j) )
    g1(12,i,j) = f1(5,i,j)+f1(6,i,j)+f1(7,i,j)+f1(8,i,j)

    g2(0,i,j) = f2(0,i,j)+f2(1,i,j)+f2(2,i,j)+f2(3,i,j)+f2(4,i,j)+f2(5,i,j)+f2(6,i,j) &
                +f2(7,i,j)+f2(8,i,j)+f2(9,i,j)+f2(10,i,j)+f2(11,i,j)+f2(12,i,j)
    g2(1,i,j) = f2(1,i,j)-f2(3,i,j)+f2(5,i,j)-f2(6,i,j)-f2(7,i,j)+f2(8,i,j)+2.0d0*f2(9,i,j)-2.0d0*f2(11,i,j)
    g2(2,i,j) = f2(2,i,j)-f2(4,i,j)+f2(5,i,j)+f2(6,i,j)-f2(7,i,j)-f2(8,i,j)+2.0d0*f2(10,i,j)-2.0d0*f2(12,i,j)
    g2(3,i,j) = f2(5,i,j)-f2(6,i,j)+f2(7,i,j)-f2(8,i,j)
    g2(4,i,j) = f2(1,i,j)+f2(2,i,j)+f2(3,i,j)+f2(4,i,j)+2.0d0*( f2(5,i,j)+f2(6,i,j)+f2(7,i,j)+f2(8,i,j) ) &
                +4.0d0*( f2(9,i,j)+f2(10,i,j)+f2(11,i,j)+f2(12,i,j) )
    g2(5,i,j) = f2(1,i,j)-f2(2,i,j)+f2(3,i,j)-f2(4,i,j)+4.0d0*( f2(9,i,j)-f2(10,i,j)+f2(11,i,j)-f2(12,i,j) )
    g2(6,i,j) = f2(1,i,j)-f2(3,i,j)+2.0d0*( f2(5,i,j)-f2(6,i,j)-f2(7,i,j)+f2(8,i,j) )+8.0d0*( f2(9,i,j)-f2(11,i,j))
    g2(7,i,j) = f2(2,i,j)-f2(4,i,j)+2.0d0*( f2(5,i,j)+f2(6,i,j)-f2(7,i,j)-f2(8,i,j) )+8.0d0*( f2(10,i,j)-f2(12,i,j) )
    g2(8,i,j) = f2(1,i,j)-f2(3,i,j)+f2(5,i,j)-f2(6,i,j)-f2(7,i,j)+f2(8,i,j)+8.0d0*( f2(9,i,j)-f2(11,i,j) )
    g2(9,i,j) = f2(2,i,j)-f2(4,i,j)+f2(5,i,j)+f2(6,i,j)-f2(7,i,j)-f2(8,i,j)+8.0d0*( f2(10,i,j)-f2(12,i,j) )
    g2(10,i,j) = f2(1,i,j)+f2(2,i,j)+f2(3,i,j)+f2(4,i,j)+2.0d0*( f2(5,i,j)+f2(6,i,j)+f2(7,i,j)+f2(8,i,j) ) &
                + 16.0d0*( f2(9,i,j)+f2(10,i,j)+f2(11,i,j)+f2(12,i,j) )
    g2(11,i,j) = f2(1,i,j)-f2(2,i,j)+f2(3,i,j)-f2(4,i,j)+16.0d0*( f2(9,i,j)-f2(10,i,j)+f2(11,i,j)-f2(12,i,j) )
    g2(12,i,j) = f2(5,i,j)+f2(6,i,j)+f2(7,i,j)+f2(8,i,j)

            s1(0) = 0.0d0
            s1(1) = 1.0d0/tau
            s1(2) = 1.0d0/tau
            s1(3) = 1.0d0/tau1
            s1(4) = 1.0d0/tau1
            s1(5) = 1.0d0/tau1
            s1(6) = 1.0d0/tau1
            s1(7) = 1.0d0/tau1
            s1(8) = 1.0d0/tau1
            s1(9) = 1.0d0/tau1
            s1(10) = 1.0d0/tau1
            s1(11) = 1.0d0/tau1
            s1(12) = 1.0d0/tau1

            s2(0) = 0.0d0
            s2(1) = 1.0d0/tau
            s2(2) = 1.0d0/tau
            s2(3) = 1.0d0/tau2
            s2(4) = 1.0d0/tau2
            s2(5) = 1.0d0/tau2
            s2(6) = 1.0d0/tau2
            s2(7) = 1.0d0/tau2
            s2(8) = 1.0d0/tau2
            s2(9) = 1.0d0/tau2
            s2(10) = 1.0d0/tau2
            s2(11) = 1.0d0/tau2
            s2(12) = 1.0d0/tau2

            do alpha=0,12
                g1(alpha,i,j) = g1(alpha,i,j)-s1(alpha)*( g1(alpha,i,j)-geq(alpha,rho1(i,j),r1,u(i,j),v(i,j)) )
                g2(alpha,i,j) = g2(alpha,i,j)-s2(alpha)*( g2(alpha,i,j)-geq(alpha,rho2(i,j),r2,u(i,j),v(i,j)) )
            enddo

    f1(0,i,j) = ( 4.0d0*g1(0,i,j)-5.0d0*g1(4,i,j)+g1(10,i,j)+4.0d0*g1(12,i,j) )/4.0d0
    f1(1,i,j) = ( 8.0d0*g1(1,i,j)+4.0d0*g1(4,i,j)+4.0d0*g1(5,i,j)-6.0d0*g1(6,i,j) &
                +4.0d0*g1(8,i,j)-g1(10,i,j)-g1(11,i,j)-6.0d0*g1(12,i,j) )/12.0d0
    f1(2,i,j) = ( 8.0d0*g1(2,i,j)+4.0d0*g1(4,i,j)-4.0d0*g1(5,i,j)-6.0d0*g1(7,i,j) &
                +4.0d0*g1(9,i,j)-g1(10,i,j)+g1(11,i,j)-6.0d0*g1(12,i,j) )/12.0d0
    f1(3,i,j) = ( -8.0d0*g1(1,i,j)+4.0d0*g1(4,i,j)+4.0d0*g1(5,i,j)+6.0d0*g1(6,i,j) &
                -4.0d0*g1(8,i,j)-g1(10,i,j)-g1(11,i,j)-6.0d0*g1(12,i,j) )/12.0d0
    f1(4,i,j) = ( -8.0d0*g1(2,i,j)+4.0d0*g1(4,i,j)-4.0d0*g1(5,i,j)+6.0d0*g1(7,i,j) &
                -4.0d0*g1(9,i,j)-g1(10,i,j)+g1(11,i,j)-6.0d0*g1(12,i,j) )/12.0d0
    f1(5,i,j) = ( g1(3,i,j)+g1(6,i,j)+g1(7,i,j)-g1(8,i,j)-g1(9,i,j)+g1(12,i,j) )/4.0d0
    f1(6,i,j) = ( -g1(3,i,j)-g1(6,i,j)+g1(7,i,j)+g1(8,i,j)-g1(9,i,j)+g1(12,i,j) )/4.0d0
    f1(7,i,j) = ( g1(3,i,j)-g1(6,i,j)-g1(7,i,j)+g1(8,i,j)+g1(9,i,j)+g1(12,i,j) )/4.0d0
    f1(8,i,j) = ( -g1(3,i,j)+g1(6,i,j)-g1(7,i,j)-g1(8,i,j)+g1(9,i,j)+g1(12,i,j) )/4.0d0
    f1(9,i,j) = ( -4.0d0*g1(1,i,j)-g1(4,i,j)-g1(5,i,j)+4.0d0*g1(8,i,j)+g1(10,i,j)+g1(11,i,j) )/48.0d0
    f1(10,i,j) = ( -4.0d0*g1(2,i,j)-g1(4,i,j)+g1(5,i,j)+4.0d0*g1(9,i,j)+g1(10,i,j)-g1(11,i,j) )/48.0d0
    f1(11,i,j) = ( 4.0d0*g1(1,i,j)-g1(4,i,j)-g1(5,i,j)-4.0d0*g1(8,i,j)+g1(10,i,j)+g1(11,i,j) )/48.0d0
    f1(12,i,j) = ( 4.0d0*g1(2,i,j)-g1(4,i,j)+g1(5,i,j)-4.0d0*g1(9,i,j)+g1(10,i,j)-g1(11,i,j) )/48.0d0

    f2(0,i,j) = ( 4.0d0*g2(0,i,j)-5.0d0*g2(4,i,j)+g2(10,i,j)+4.0d0*g2(12,i,j) )/4.0d0
    f2(1,i,j) = ( 8.0d0*g2(1,i,j)+4.0d0*g2(4,i,j)+4.0d0*g2(5,i,j)-6.0d0*g2(6,i,j) &
                +4.0d0*g2(8,i,j)-g2(10,i,j)-g2(11,i,j)-6.0d0*g2(12,i,j) )/12.0d0
    f2(2,i,j) = ( 8.0d0*g2(2,i,j)+4.0d0*g2(4,i,j)-4.0d0*g2(5,i,j)-6.0d0*g2(7,i,j) &
                +4.0d0*g2(9,i,j)-g2(10,i,j)+g2(11,i,j)-6.0d0*g2(12,i,j) )/12.0d0
    f2(3,i,j) = ( -8.0d0*g2(1,i,j)+4.0d0*g2(4,i,j)+4.0d0*g2(5,i,j)+6.0d0*g2(6,i,j) &
                -4.0d0*g2(8,i,j)-g2(10,i,j)-g2(11,i,j)-6.0d0*g2(12,i,j) )/12.0d0
    f2(4,i,j) = ( -8.0d0*g2(2,i,j)+4.0d0*g2(4,i,j)-4.0d0*g2(5,i,j)+6.0d0*g2(7,i,j) &
                -4.0d0*g2(9,i,j)-g2(10,i,j)+g2(11,i,j)-6.0d0*g2(12,i,j) )/12.0d0
    f2(5,i,j) = ( g2(3,i,j)+g2(6,i,j)+g2(7,i,j)-g2(8,i,j)-g2(9,i,j)+g2(12,i,j) )/4.0d0
    f2(6,i,j) = ( -g2(3,i,j)-g2(6,i,j)+g2(7,i,j)+g2(8,i,j)-g2(9,i,j)+g2(12,i,j) )/4.0d0
    f2(7,i,j) = ( g2(3,i,j)-g2(6,i,j)-g2(7,i,j)+g2(8,i,j)+g2(9,i,j)+g2(12,i,j) )/4.0d0
    f2(8,i,j) = ( -g2(3,i,j)+g2(6,i,j)-g2(7,i,j)-g2(8,i,j)+g2(9,i,j)+g2(12,i,j) )/4.0d0
    f2(9,i,j) = ( -4.0d0*g2(1,i,j)-g2(4,i,j)-g2(5,i,j)+4.0d0*g2(8,i,j)+g2(10,i,j)+g2(11,i,j) )/48.0d0
    f2(10,i,j) = ( -4.0d0*g2(2,i,j)-g2(4,i,j)+g2(5,i,j)+4.0d0*g2(9,i,j)+g2(10,i,j)-g2(11,i,j) )/48.0d0
    f2(11,i,j) = ( 4.0d0*g2(1,i,j)-g2(4,i,j)-g2(5,i,j)-4.0d0*g2(8,i,j)+g2(10,i,j)+g2(11,i,j) )/48.0d0
    f2(12,i,j) = ( 4.0d0*g2(2,i,j)-g2(4,i,j)+g2(5,i,j)-4.0d0*g2(9,i,j)+g2(10,i,j)-g2(11,i,j) )/48.0d0

        enddo
    enddo

    return
    end subroutine collision


    function geq(alpha,rho,rk,u,v)
    implicit none
    real(8) :: geq
    integer :: alpha
    real(8) :: rk
    real(8) :: u, v, rho

    select case (alpha)
        case (0)
            geq = rho
        case (1)
            geq = rho*u
        case (2)
            geq = rho*v
        case (3)
            geq = rho*u*v
        case (4)
            geq = rho*( 8.0d0/(5.0d0+3.0d0*rk)+u*u+v*v)
        case (5)
            geq = rho*(u*u-v*v)
        case (6)
            geq = 16.0d0*rho*u/(5.0d0+3.0d0*rk)
        case (7)
            geq = 16.0d0*rho*v/(5.0d0+3.0d0*rk)
        case (8)
            geq = 12.0d0*rho*u/(5.0d0+3.0d0*rk)
        case (9)
            geq = 12.0d0*rho*v/(5.0d0+3.0d0*rk)
        case (10)
            geq = rho/2.0d0*( 24.0d0/(5.0d0+3.0d0*rk)+5.0d0*(u*u+v*v) )
        case (11)
            geq = 3.0d0*rho*( u*u-v*v)
        case (12)
            geq = rho/4.0d0*( 8.0d0/(5.0d0+3.0d0*rk)+u*u+v*v )
    end select

    return
    end function geq

!!! streaming step
    subroutine streaming()
    use alldata
    implicit none
    integer :: i, j
!    integer :: id, jd
!    integer :: alpha

        do i=1,nx
            do j=1,ny
                f1(0,i,j) = f1(0,i,j)
            enddo
        enddo
        do i=nx,2,-1
            do j=1,ny
                f1(1,i,j) = f1(1,i-1,j)
            enddo
        enddo
        do i=nx,3,-1
            do j=1,ny
                f1(9,i,j) = f1(9,i-2,j)
            enddo
        enddo
        do i=1,nx
            do j=ny,2,-1
                f1(2,i,j) = f1(2,i,j-1)
            enddo
        enddo
        do i=1,nx
            do j=ny,3,-1
                f1(10,i,j) = f1(10,i,j-2)
            enddo
        enddo
        do i=1,nx-1
            do j=1,ny
                f1(3,i,j) = f1(3,i+1,j)
            enddo
        enddo
        do i=1,nx-2
            do j=1,ny
                f1(11,i,j) = f1(11,i+2,j)
            enddo
        enddo
        do i=1,nx
            do j=1,ny-1
                f1(4,i,j) = f1(4,i,j+1)
            enddo
        enddo
        do i=1,nx
            do j=1,ny-2
                f1(12,i,j) = f1(12,i,j+2)
            enddo
        enddo
        do i=nx,2,-1
            do j=ny,2,-1
                f1(5,i,j) = f1(5,i-1,j-1)
            enddo
        enddo
        do i=1,nx-1
            do j=ny,2,-1
                f1(6,i,j) = f1(6,i+1,j-1)
            enddo
        enddo
        do i=1,nx-1
            do j=1,ny-1
                f1(7,i,j) = f1(7,i+1,j+1)
            enddo
        enddo
        do i=nx,2,-1
            do j=1,ny-1
                f1(8,i,j) = f1(8,i-1,j+1)
            enddo
        enddo

        !Periodic boundary conditon
        do j=1,ny
            f1(1,1,j) = f1(1,nx,j)
            f1(5,1,j) = f1(5,nx,j)
            f1(8,1,j) = f1(8,nx,j)
            f1(9,2,j) = f1(9,nx,j)

            f1(3,nx,j) = f1(3,1,j)
            f1(6,nx,j) = f1(6,1,j)
            f1(7,nx,j) = f1(7,1,j)
            f1(11,nx,j) = f1(11,2,j)
        enddo

        do i=1,nx
            f1(4,i,ny) = f1(4,i,1)
            f1(12,i,ny) = f1(12,i,2)

            f1(2,i,1) = f1(2,i,ny)
            f1(10,i,1) = f1(10,i,ny-1)
        enddo

!----------------------------------------

        do i=1,nx
            do j=1,ny
                f2(0,i,j) = f2(0,i,j)
            enddo
        enddo
        do i=nx,2,-1
            do j=1,ny
                f2(1,i,j) = f2(1,i-1,j)
            enddo
        enddo
        do i=nx,3,-1
            do j=1,ny
                f2(9,i,j) = f2(9,i-2,j)
            enddo
        enddo
        do i=1,nx
            do j=ny,2,-1
                f2(2,i,j) = f2(2,i,j-1)
            enddo
        enddo
        do i=1,nx
            do j=ny,3,-1
                f2(10,i,j) = f2(10,i,j-2)
            enddo
        enddo
        do i=1,nx-1
            do j=1,ny
                f2(3,i,j) = f2(3,i+1,j)
            enddo
        enddo
        do i=1,nx-2
            do j=1,ny
                f2(11,i,j) = f2(11,i+2,j)
            enddo
        enddo
        do i=1,nx
            do j=1,ny-1
                f2(4,i,j) = f2(4,i,j+1)
            enddo
        enddo
        do i=1,nx
            do j=1,ny-2
                f2(12,i,j) = f2(12,i,j+2)
            enddo
        enddo
        do i=nx,2,-1
            do j=ny,2,-1
                f2(5,i,j) = f2(5,i-1,j-1)
            enddo
        enddo
        do i=1,nx-1
            do j=ny,2,-1
                f2(6,i,j) = f2(6,i+1,j-1)
            enddo
        enddo
        do i=1,nx-1
            do j=1,ny-1
                f2(7,i,j) = f2(7,i+1,j+1)
            enddo
        enddo
        do i=nx,2,-1
            do j=1,ny-1
                f2(8,i,j) = f2(8,i-1,j+1)
            enddo
        enddo

        !Periodic boundary conditon
        do j=1,ny
            f2(1,1,j) = f2(1,nx,j)
            f2(5,1,j) = f2(5,nx,j)
            f2(8,1,j) = f2(8,nx,j)
            f2(9,2,j) = f2(9,nx,j)

            f2(3,nx,j) = f2(3,1,j)
            f2(6,nx,j) = f2(6,1,j)
            f2(7,nx,j) = f2(7,1,j)
            f2(11,nx,j) = f2(11,2,j)
        enddo

        do i=1,nx
            f2(4,i,ny) = f2(4,i,1)
            f2(12,i,ny) = f2(12,i,2)

            f2(2,i,1) = f2(2,i,ny)
            f2(10,i,1) = f2(10,i,ny-1)
        enddo

!!! Warning: This section may not be right!
!!!    do i=3,nx-2
!!!        do j=3,ny-2
!!!            do alpha=0,12
!!!                id = MOD( INT(i-ex(alpha)+nx),nx )
!!!                jd = MOD( INT(j-ey(alpha)+ny),ny )
!!!                f1(alpha,i,j) = f1(alpha,id,jd)
!!!                f2(alpha,i,j) = f2(alpha,id,jd)
!!!            enddo
!!!        enddo
!!!    enddo
!!! Warning: This section may not be right!

    return
    end subroutine streaming


!!! check convergence
    subroutine check()
    use alldata
    implicit none
    integer :: i, j

    itc = itc+1
    error = 0.0d0
    if(itc.EQ.1) error = 10.0d0
    if(itc.EQ.2) error = 10.0d0
    if(itc.EQ.3) error = 10.0d0

    if(itc.GT.3) then
        do i=1,nx
            do j=1,ny
                error  = error+SQRT((u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))) &
                                /SQRT((u(i,j)+0.00001)*(u(i,j)+0.00001)+(v(i,j)+0.00001)*(v(i,j)+0.00001))
            enddo
        enddo
    endif

    up = u
    vp = v

    if(MOD(itc,50).EQ.0)  write(*,*) itc,' ',error

!!!        open(unit=01,file='error.dat',status='unknown',position='append')
!!!        if (MOD(itc,2000).EQ.0) then
!!!            write(01,*) itc,' ',error
!!!        endif
!!!        close(01)

    return
    end subroutine check

!!! output data file
    subroutine output()
    use alldata
    implicit none
    integer :: i, j
    character(len=100) :: filename
    write(filename,*) itc
    filename = adjustl(filename)

    open(unit=02,file='sineWave-'//trim(filename)//'.plt',status='unknown')

    write(02,*) 'TITLE="Decay of the sine-wave density profile"'
    write(02,*) 'VARIABLES="RHO1" "RHO2"'
    write(02,101) nx, ny
    do j=1,ny
        do i = 1,nx
            write(02,100) rho1(i,j), rho2(i,j)
        enddo
    enddo

100 format(2x,10(e12.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'DATAPACKING=POINT')

    return
    end subroutine output
