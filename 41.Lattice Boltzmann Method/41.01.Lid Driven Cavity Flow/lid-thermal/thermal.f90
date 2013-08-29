
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
        call streaming_f()

!!! boundary condition
        call bounceback_f()

!!! collision step
        call collision_f()

!!! streaming step
        call streaming_T()

!!! boundary condition
        call bounceback_T()

!!! collision step
        call collision_T()

!!! check convergence
        call check(up,vp,error)

!!! output preliminary results
        if(MOD(itc,500).EQ.0) then
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
    write(*,*) 'Ra =',Ra
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
    real(8) :: un(0:8)

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
    temperature = 0.0d0
    do j=1,ny
        temperature(nx,j) = 1.0d0
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
    subroutine streaming_f()
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j

    do i=1,nx
        do j=1,ny
            f(0,i,j) = f(0,i,j)
        enddo
    enddo

    do i=nx,2,-1
        do j=1,ny
            f(1,i,j) = f(1,i-1,j)
        enddo
    enddo

    do i=1,nx
        do j=ny,2,-1
            f(2,i,j) = f(2,i,j-1)
        enddo
    enddo

    do i=1,nx-1
        do j=1,ny
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
        do j=ny,2,-1
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
    end subroutine streaming_f

!!! collision step
    subroutine collision_f()
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j
    integer :: alpha
    real(8) :: us2
    real(8) :: un(0:8)
    real(8) :: feq(0:8,nx,ny)

    do i=1,nx
        do j=1,ny

            rho(i,j) = 0.0d0
            do alpha=0,8
                rho(i,j) = rho(i,j)+f(alpha,i,j)
            enddo

            do alpha=0,8
                force(alpha,i,j) = 0.0d0
            enddo
            force(2,i,j) = 0.5d0*1.0d0/3.0d0*Ma*Ma*temperature(i,j)
            force(4,i,j) = -0.5d0*1.0d0/3.0d0*Ma*Ma*temperature(i,j)

            !data ex/0.0d0,1.0d0,0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0/
            !data ey/0.0d0,0.0d0,1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0/
            u(i,j) = (f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j))/rho(i,j)
            v(i,j) = (f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j))/rho(i,j)
            us2 = u(i,j)*u(i,j)+v(i,j)*v(i,j)

            do alpha=0,8
                un(alpha) = u(i,j)*ex(alpha) + v(i,j)*ey(alpha)
                feq(alpha,i,j) = omega(alpha)*rho(i,j) &
                            *(1.0d0+un(alpha)/cs2+un(alpha)*un(alpha)/(2.0d0*cs2*cs2)-us2/(2.0d0*cs2))
                f(alpha,i,j) = f(alpha,i,j)-1.0d0/tau_f*(f(alpha,i,j)-feq(alpha,i,j)) &
                            +dt*force(alpha,i,j)
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
    end subroutine collision_f

!!! boundary condition
    subroutine bounceback_f()
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

    do i=1,nx
        !Bottom side
        f(2,i,1) = f(4,i,1)
        f(5,i,1) = f(7,i,1)
        f(6,i,1) = f(8,i,1)

        !Top side
        f(4,i,nx) = f(2,i,nx)
        f(7,i,nx) = f(5,i,nx)
        f(8,i,nx) = f(6,i,nx)

    enddo

!    !Left-Bottom corner
!    f(1,1,1) = f(3,1,1)
!    f(2,1,1) = f(4,1,1)
!    f(5,1,1) = f(7,1,j)
!
!    !Right-Bottom corner
!    f(3,nx,1) = f(1,nx,1)
!    f(6,nx,1) = f(8,nx,1)
!    f(2,nx,1) = f(4,nx,1)

    return
    end subroutine bounceback_f



!!! streaming step
    subroutine streaming_T()
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j

!    do i=1,nx
!        do j=1,ny-1
!            T(0,i,j) = T(0,i,j)
!        enddo
!    enddo

    do i=nx,2,-1
        do j=1,ny
            T(1,i,j) = T(1,i-1,j)
        enddo
    enddo

    do i=1,nx
        do j=ny,2,-1
            T(2,i,j) = T(2,i,j-1)
        enddo
    enddo

    do i=2,nx-1
        do j=1,ny
            T(3,i,j) = T(3,i+1,j)
        enddo
    enddo

    do i=1,nx
        do j=1,ny-1
            T(4,i,j) = T(4,i,j+1)
        enddo
    enddo

!    do i=nx,2,-1
!        do j=ny-1,2,-1
!            T(5,i,j) = T(5,i-1,j-1)
!        enddo
!    enddo
!
!    do i=1,nx-1
!        do j=ny-1,2,-1
!            T(6,i,j) = T(6,i+1,j-1)
!        enddo
!    enddo
!
!    do i=1,nx-1
!        do j=1,ny-1
!            T(7,i,j) = T(7,i+1,j+1)
!        enddo
!    enddo
!
!    do i=nx,2,-1
!        do j=1,ny-1
!            T(8,i,j) = T(8,i-1,j+1)
!        enddo
!    enddo

    return
    end subroutine streaming_T

!!! collision step
    subroutine collision_T()
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j
    integer :: alpha
    real(8) :: un(0:8)
    real(8) :: Teq(0:4,nx,ny)

    do i=2,nx-1
        do j=1,ny

            temperature(i,j) = 0.0d0
            do alpha = 1,4
                temperature(i,j) = temperature(i,j)+T(alpha,i,j)
            enddo

            do alpha=1,4
                un(alpha) = u(i,j)*ex(alpha) + v(i,j)*ey(alpha)
!                Teq(alpha,i,j) = omega(alpha)*temperature(i,j) &
!                            *(1.0d0+un(alpha)/cs2+4.5d0*un(alpha)*un(alpha)/(cs2*cs2)-1.5d0*us2/cs2)
                Teq(alpha,i,j) = temperature(i,j)/4.0d0*(1.0d0+2.0d0*un(alpha)/cs2)
                T(alpha,i,j) = T(alpha,i,j)-1.0d0/tau_T*(T(alpha,i,j)-Teq(alpha,i,j))
            enddo

        enddo
    enddo

    do j=1,ny
        temperature(1,j) = 1.0d0
        temperature(nx,j) = 0.0d0
    enddo


    return
    end subroutine collision_T

!!! boundary condition
    subroutine bounceback_T()
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j

    do j=2,ny-1
        !Left side
        T(1,1,j) = T(3,1,j)
!        T(5,1,j) = T(7,1,j)
!        T(8,1,j) = T(6,1,j)

        !Right side
        T(3,nx,j) = T(1,nx,j)
!        T(6,nx,j) = T(8,nx,j)
!        T(7,nx,j) = T(5,nx,j)
    enddo

    do i=1,nx
        !Bottom side
        T(2,i,1) = T(4,i,1)
!        T(5,i,1) = T(7,i,1)
!        T(6,i,1) = T(8,i,1)

        !Top side
        T(4,i,ny) = T(2,i,ny)
!        T(7,i,ny) = T(5,i,ny)
!        T(8,i,ny) = T(6,i,ny)
    enddo

!    !Left-Bottom corner
!    T(1,1,1) = T(3,1,1)
!    T(2,1,1) = T(4,1,1)
!    T(5,1,1) = T(7,1,j)
!
!    !Right-Bottom corner
!    T(3,nx,1) = T(1,nx,1)
!    T(6,nx,1) = T(8,nx,1)
!    T(2,nx,1) = T(4,nx,1)

    return
    end subroutine bounceback_T
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
        do j=1,ny
            p(i,j) = rho(i,j)*cs2
        enddo
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

        do j=1,ny
            psi(1,j) = 0.0d0
            psi(nx,j) = 0.0d0
        enddo
        do i=1,nx
            psi(i,1) = 0.0d0
            psi(i,ny) = 0.0d0
        enddo

    do i=3,nx-2
            psi(i,3) = up(i,2)*2.0d0*dy+psi(i,1)
            psi(i,2) = 0.25d0*psi(i,3)
        do j=3,ny-3
            psi(i,j+1) = up(i,j)*2.0d0*dy+psi(i,j-1)
            !psi(i+1,j) = -v(i,j)*2.0d0*dx+psi(i-1,j) ! Alternative and equivalent psi formulae
        enddo
        psi(i,ny-1) = 0.25d0*psi(i,ny-2)
    enddo

    do j=2,ny-1
        psi(2,j) = 0.25d0*psi(3,j)
        psi(nx-1,j) = 0.25d0*psi(nx-2,j)
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

    open(unit=02,file='rb_'//trim(filename)//'.plt',status='unknown')
    write(02,*) 'TITLE="Lid Driven Cavity(LBM)"'
    write(02,*) 'VARIABLES=x,y,u,v,psi,temperature'
    write(02,101) nx, ny
    do j=1,ny
        do i = 1,nx
            write(02,100) X(i), Y(j), up(i,j), vp(i,j), psi(i,j), temperature(i,j)
        enddo
    enddo

100 format(2x,10(e12.6,' '))
101 format('ZONE',1x,'i=',1x,i5,2x,'j=',1x,i5,1x,'F=POINT')

    close(02)

    return
    end subroutine output
