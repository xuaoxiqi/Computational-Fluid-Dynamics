
!!!    This program sloves Lid Driven Cavity Flow problem using Lattice Boltzmann Method
!!!    Lattice Boltzmann Equation with BGK approximation
!!!    Copyright (C) 2013  Ao Xu
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

!!!                  Moving Wall
!!!               |---------------|
!!!               |               |
!!!               |               |
!!!    Stationary |               | Stationary
!!!       Wall    |               |    Wall
!!!               |               |
!!!               |               |
!!!               |---------------|
!!!                Stationary Wall

!!!     D2Q9 Lattice Vector Properties:
!!!              6   2   5
!!!                \ | /
!!!              3 - 0 - 1
!!!                / | \
!!!              7   4   8


    program main
    implicit none
    include "para.h"
    include "common.h"
    real(8) :: error = 100.0d0
    real(8) :: up(nx,ny), vp(nx,ny)

!!! set up initial flow field
    call initial()

    do while((error.GT.eps).AND.(itc.LT.itc_max))

!!! streaming step
        call streaming()

!!! boundary condition
        call bounceback()

!!! collision step
        call collision()

!!! check convergence
        call check(up,vp,error)

!!! output preliminary results
        if(MOD(itc,2000).EQ.0) then
            call calp()
            call calpsi(up,vp)
            call output(up,vp)
        endif

    enddo

!!! compute pressure field
    call calp()

!!! compute streamfunction
    call calpsi(up,vp)

!!! output data file
    call output(up,vp)

    write(*,*)
    write(*,*) '************************************************************'
    write(*,*) 'This program sloves Lid Driven Cavity Flow problem using Lattice Boltzmann Method'
    write(*,*) 'Lattice Boltzmann Equation with BGK approximation'
    write(*,*) 'Consider D2Q9 Particle Discrete Velocity model'
    write(*,*) 'nx =',nx,',       ny =',ny
    write(*,*) 'Re =',Re
    write(*,*) 'eps =',eps
    write(*,*) 'itc =',itc
    write(*,*) '************************************************************'
    write(*,*)

    stop
    end program main

!!! set up initial flow field
    subroutine initial()
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j
    integer :: alpha
    real(8) :: us2

    do i=1,nx
        X(i) = (i-1)*dx
    enddo
    do j=1,ny
        Y(j) = (j-1)*dy
    enddo

    omega(0) = 4.0d0/9.0d0
    do alpha=1,4
        omega(alpha) = 1.0d0/9.0d0
    enddo
    do alpha=5,8
        omega(alpha) = 1.0d0/36.0d0
    enddo

    u = 0.0d0
    v = 0.0d0
    rho = 1.0d0
    psi = 0.0d0
    do i=1,nx
        u(i,ny) = U_ref
    enddo

    itc = 0

    do i=1,nx
        do j=1,ny
            us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha=0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                f(alpha,i,j) = omega(alpha)*(1.0d0+un(alpha)/cs2+un(alpha)*un(alpha)/(2.0d0*cs2*cs2)-us2/(2.0d0*cs2))
            enddo
        enddo
    enddo

    return
    end subroutine initial

!!! streaming step
    subroutine streaming()
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j

    do i=1,nx
        do j=1,ny-1
            f(0,i,j) = f(0,i,j)
        enddo
    enddo

    do i=nx,2,-1
        do j=1,ny-1
            f(1,i,j) = f(1,i-1,j)
        enddo
    enddo

    do i=1,nx
        do j=ny-1,2,-1
            f(2,i,j) = f(2,i,j-1)
        enddo
    enddo

    do i=1,nx-1
        do j=1,ny-1
            f(3,i,j) = f(3,i+1,j)
        enddo
    enddo

    do i=1,nx
        do j=1,ny-1
            f(4,i,j) = f(4,i,j+1)
        enddo
    enddo

    do i=nx,2,-1
        do j=ny-1,2,-1
            f(5,i,j) = f(5,i-1,j-1)
        enddo
    enddo

    do i=1,nx-1
        do j=ny-1,2,-1
            f(6,i,j) = f(6,i+1,j-1)
        enddo
    enddo

    do i=1,nx-1
        do j=1,ny-1
            f(7,i,j) = f(7,i+1,j+1)
        enddo
    enddo

    do i=nx,2,-1
        do j=1,ny-1
            f(8,i,j) = f(8,i-1,j+1)
        enddo
    enddo

    return
    end subroutine streaming

!!! collision step
    subroutine collision()
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j
    integer :: alpha
    real(8) :: us2
    real(8) :: feq(0:8,nx,ny)

    do i=1,nx
        do j=1,ny-1

            rho(i,j) = 0.0d0
            do alpha=0,8
                rho(i,j) = rho(i,j)+f(alpha,i,j)
            enddo

            !data ex/0.0d0,1.0d0,0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0/
            !data ey/0.0d0,0.0d0,1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0/
            u(i,j) = (f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j))/rho(i,j)
            v(i,j) = (f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j))/rho(i,j)
            us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)

            do alpha=0,8
                un(alpha) = u(i,j)*ex(alpha) + v(i,j)*ey(alpha)
                feq(alpha,i,j) = omega(alpha)*rho(i,j) &
                            *(1.0d0+un(alpha)/cs2+un(alpha)*un(alpha)/(2.0d0*cs2*cs2)-us2/(2.0d0*cs2))
                f(alpha,i,j) = f(alpha,i,j)-1.0d0/tau*(f(alpha,i,j)-feq(alpha,i,j))
            enddo

        enddo
    enddo

    !Left bottom corner
    f(6,1,1) = feq(6,1,1)
    f(8,1,1) = feq(8,1,1)

    !Right bottom corner
    f(5,nx,1) = feq(5,nx,1)
    f(7,nx,1) = feq(7,nx,1)

    return
    end subroutine collision

!!! boundary condition
    subroutine bounceback()
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j

    do j=2,ny-1
        !Left side
        f(1,1,j) = f(3,1,j)
        f(5,1,j) = f(7,1,j)
        f(8,1,j) = f(6,1,j)

        !Right side
        f(3,nx,j) = f(1,nx,j)
        f(6,nx,j) = f(8,nx,j)
        f(7,nx,j) = f(5,nx,j)
    enddo

    do i=2,nx-1
        !Bottom side
        f(2,i,1) = f(4,i,1)
        f(5,i,1) = f(7,i,1)
        f(6,i,1) = f(8,i,1)
    enddo

    !Left-Bottom corner
    f(1,1,1) = f(3,1,1)
    f(2,1,1) = f(4,1,1)
    f(5,1,1) = f(7,1,j)

    !Right-Bottom corner
    f(3,nx,1) = f(1,nx,1)
    f(6,nx,1) = f(8,nx,1)
    f(2,nx,1) = f(4,nx,1)

    return
    end subroutine bounceback

!!! check convergence
    subroutine check(up,vp,error)
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j
    real(8) :: error
    real(8) :: up(nx,ny), vp(nx,ny)

    itc = itc+1
    error = 0.0d0
    if(itc.EQ.1) error = 10.0d0
    if(itc.EQ.2) error = 10.0d0
    if(itc.EQ.3) error = 10.0d0

    if(itc.GT.3) then
        do i=1,nx
            do j=1,ny-1
                error  = error+SQRT((u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))) &
                        /SQRT((u(i,j)+0.00001)*(u(i,j)+0.00001)+(v(i,j)+0.00001)*(v(i,j)+0.00001))
            enddo
        enddo
    endif

    up = u
    vp = v

    if(MOD(itc,50).EQ.0) write(*,*) itc,' ',error

!!!        open(unit=01,file='error.dat',status='unknown',position='append')
!!!        if (MOD(itc,2000).EQ.0) then
!!!            write(01,*) itc,' ',error
!!!        endif
!!!        close(01)

    return
    end subroutine check

!!! compute pressure field
    subroutine calp()
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j

    do i=1,nx
        do j=1,ny-1
            p(i,j) = rho(i,j)*cs2
        enddo
    enddo

    do i=1,nx
        p(i,ny) = cs2
    enddo

    return
    end subroutine calp

!!! compute Streamfunction
    subroutine calpsi(up,vp)
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j
    real(8) :: up(nx,ny), vp(nx,ny)

!        do j=1,ny
!            psi(1,j) = 0.0d0
!            psi(nx,j) = 0.0d0
!        enddo
!        do i=1,nx
!            psi(i,1) = 0.0d0
!            psi(i,ny) = 0.0d0
!        enddo

    do i=3,nx-2
        do j=2,ny-3
            psi(i,j+1) = up(i,j)*2.0d0*dy+psi(i,j-1)
            !psi(i+1,j) = -v(i-1,j)*2.0d0*dx+psi(i-1,j) ! Alternative and equivalent psi formulae
        enddo
    enddo

    do j=2,ny-1
        psi(2,j) = 0.25d0*psi(3,j)
        psi(nx-1,j) = 0.25d0*psi(nx-2,j)
    enddo
    do i=2,nx-1
        psi(i,2) = 0.25d0*psi(i,3)
        psi(i,ny-1) = 0.25d0*(psi(i,ny-2)-0.2d0*dy)
    enddo

    return
    end subroutine calpsi

!!! output data file
    subroutine output(up,vp)
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j
    real(8) :: up(nx,ny), vp(nx,ny)
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)

    open(unit=02,file='cavity_'//trim(filename)//'.plt',status='unknown')
    write(02,*) 'TITLE="Lid Driven Cavity(LBM)"'
    write(02,*) 'VARIABLES=x,y,u,v,psi,p'
    write(02,101) nx, ny
    do j=1,ny
        do i = 1,nx
            write(02,100) X(i), Y(j), up(i,j), vp(i,j), psi(i,j), p(i,j)
        enddo
    enddo

100 format(2x,10(e12.6,' '))
101 format('ZONE',1x,'i=',1x,i5,2x,'j=',1x,i5,1x,'F=POINT')

    close(02)

    return
    end subroutine output
