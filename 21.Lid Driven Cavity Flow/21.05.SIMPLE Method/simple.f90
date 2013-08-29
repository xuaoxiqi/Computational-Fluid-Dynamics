
!!!    This program sloves Lid Driven Cavity Flow problem using SIMPLE Method
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


        program main
        implicit none
        include "para.h"
        include "common.h"
        real(8) :: un(nx,ny+1),vn(nx+1,ny),pn(nx+1,ny+1),uc(nx,ny),vc(nx,ny),pc(nx,ny)
        real(8) :: dp(nx+1,ny+1), ur(nx,ny+1), vr(nx+1,ny)
        real(8) :: D(nx+1,ny+1)
        real(8) :: error

!!! input initial data
        itc = 0
        error=100.00d0


!!! set up initial flow field
        call initial(dp)

        do while((error.GT.eps).AND.(itc.LT.itc_max))

!!! Use guessed pressure values to solve for velocity from momentum equations
            call solMom(nx,ny,dx,dy,dt,Re,u,v,pGuess,ur,vr)

!!! Use continuity equation to construct a pressure correction dp
            call caldp(nx,ny,dx,dy,dt,Re,ur,vr,dp)

!!! Update velocity using velocity correction du, dv
            call update_uvp(nx,ny,dx,dy,dt,ur,vr,pGuess,dp,un,vn,pn)

!!! check convergence
            call check(nx,ny,dx,dy,dt,error,ur,vr,un,vn,pn,u,v,pGuess,itc)

!!! output preliminary results
            if (MOD(itc,20000).EQ.0) then
                call caluvp(uc,vc,pc)
                call calpsi(uc,vc)
                call output(uc,vc)
            endif

        enddo

!!! compute velocity components u, v and pressure p
        call caluvp(uc,vc,pc)

!!! compute Streamfunction
        call calpsi(uc,vc)

!!! output data file
        call output(uc,vc)

        write(*,*)
        write(*,*) '************************************************************'
        write(*,*) 'This program sloves Lid Driven Cavity Flow problem'
        write(*,*) 'using SIMPLE Method'
        write(*,*) 'nx =',nx,',       ny =',ny
        write(*,*) 'Re =',Re
        write(*,*) 'dt =',dt
        write(*,*) 'eps =',eps
        write(*,*) 'itc =',itc
        write(*,*) 'Developing time=',dt*itc,'s'
        write(*,*) '************************************************************'
        write(*,*)

        stop
        end program main


!!! set up initial flow field
        subroutine initial(dp)
        implicit none
        include "para.h"
        include "common.h"
        integer :: i, j
        real(8) :: dp(nx+1,ny+1)

        do i=1,nx
            X(i) = (i-1)*dx
        enddo
        do j=1,ny
            Y(j) = (j-1)*dy
        enddo
        psi = 0.0d0

        pGuess = 1.0d0
        uGuess = 0.0d0
        vGuess = 0.0d0

        dp = 0.0d0

        do i=1,nx
            uGuess(i,ny+1) = 4.0d0/3.0d0
            uGuess(i,ny) = 2.0d0/3.0d0
        enddo

        return
        end subroutine initial


!!! Use guessed pressure values to solve for velocity from momentum equations
        subroutine solmom(nx,ny,dx,dy,dt,Re,u,v,pGuess,ur,vr)
        implicit none
        integer :: nx, ny, i, j
        real(8) :: u(nx,ny+1),v(nx+1,ny),pGuess(nx+1,ny+1),ur(nx,ny+1),vr(nx+1,ny)
        real(8) :: Re, dx, dy, dt

        do i=2,nx-1
            do j=2,ny
    ur(i,j) = u(i,j) - dt*(  (u(i+1,j)*u(i+1,j)-u(i-1,j)*u(i-1,j))/2.0d0/dx &
    +0.25d0*( (u(i,j)+u(i,j+1))*(v(i,j)+v(i-1,j))-(u(i,j)+u(i,j-1))*(v(i-1,j-1)+v(i,j-1)) )/dy  )&
    - dt/dx*(pGuess(i+1,j)-pGuess(i,j)) &
    + dt*1.0d0/Re*( (u(i+1,j)-2.0d0*u(i,j)+u(i-1,j))/dx/dx +(u(i,j+1)-2.0d0*u(i,j)+u(i,j-1))/dy/dy )
            enddo
        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=2,nx
            do j=2,ny-1
    vr(i,j) = v(i,j) - dt* ( 0.25d0*( (u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j))-(u(i-1,j)+u(i-1,j+1))*(v(i,j)+v(i-1,j)) )/dx &
    +(v(i,j+1)*v(i,j+1)-v(i,j-1)*v(i,j-1))/2.0d0/dy ) &
    - dt/dy*(pGuess(i,j+1)-pGuess(i,j)) &
    + dt*1.0d0/Re*( (v(i+1,j)-2.0d0*v(i,j)+v(i-1,j))/dx/dx+(v(i,j+1)-2.0d0*v(i,j)+v(i,j-1))/dy/dy )
            enddo
        enddo

        return
        end subroutine solmom


!!! Use continuity equation to construct a pressure correction dp
        subroutine caldp(nx,ny,dx,dy,dt,Re,ur,vr,dp)
        implicit none
        integer :: nx, ny, i, j
        real(8) :: Re
        real(8) :: dx, dy, dt
        real(8) :: alpha
        real(8) :: ur(nx,ny+1), vr(nx+1,ny), dp(nx+1,ny+1)
        real(8) :: R(nx,ny)

        dp = 0.0d0
        do i=2,nx
            do j=2,ny
                dp(i,j) = ( dt/dx/dx*(dp(i+1,j)+dp(i-1,j))+dt/dy/dy*(dp(i,j+1)+dp(i,j-1)) &
                            -(ur(i,j)-ur(i-1,j))/dx-(vr(i,j)-vr(i,j-1))/dy )/2.0d0/(dt/dx/dx+dt/dy/dy)
            enddo
        enddo

        return
        end subroutine caldp


!!! Update velocity using velocity correction du, dv
        subroutine update_uvp(nx,ny,dx,dy,dt,ur,vr,pGuess,dp,un,vn,pn)
        implicit none
        integer :: nx, ny, i, j
        real(8) :: dx, dy, dt
        real(8) :: du(nx,ny+1), dv(nx+1,ny), dp(nx+1,ny+1)
        real(8) :: ur(nx,ny+1), vr(nx+1,ny), pGuess(nx+1,ny+1)
        real(8) :: un(nx,ny+1), vn(nx+1,ny), pn(nx+1,ny+1)
        real(8) :: alpha

        do i=2,nx-1
            do j=2,ny
                du(i,j) = -dt/dx*(dp(i+1,j)-dp(i,j))
                un(i,j) = ur(i,j)+du(i,j)
            enddo
        enddo

        do i=2,nx
            do j=2,ny-1
                dv(i,j) = -dt/dy*(dp(i,j+1)-dp(i,j))
                vn(i,j) = vr(i,j)+dv(i,j)
            enddo
        enddo

        alpha = 0.8d0
        do i=2,nx
            do j=2,ny
                pn(i,j) = pGuess(i,j)+alpha*dp(i,j)
            enddo
        enddo

        !!! boundary condition for velocity

        do j=2,ny
            un(1,j) = 0.0d0
            un(nx,j) = 0.0d0
        enddo

        do i=1,nx
            un(i,1) = -un(i,2)
            un(i,ny+1) = 2.0d0-un(i,ny)
        enddo

        do j=2,ny-1
            vn(1,j) = -vn(2,j)
            vn(ny+1,j) = -vn(ny,j)
        enddo

        do i=1,nx+1
            vn(i,1) = 0.0d0
            vn(i,ny) = 0.0d0
        enddo

        !!! boundary condition for pressure
        do i=2,nx
            pn(i,1) = pn(i,2)
            pn(i,ny+1) = pn(i,ny)
        enddo
        do j=1,ny+1
            pn(1,j) = pn(2,j)
            pn(nx+1,j) = pn(nx,j)
        enddo

        return
        end subroutine update_uvp


!!! check convergence
        subroutine check(nx,ny,dx,dy,dt,error,ur,vr,un,vn,pn,u,v,pGuess,itc)
        implicit none
        integer :: nx, ny, i, j, itc
        real(8) :: dt, error
        real(8) :: dx, dy
        real(8) :: ur(nx,ny+1), vr(nx+1,ny)
        real(8) :: u(nx,ny+1), v(nx+1,ny), pGuess(nx+1,ny+1), un(nx,ny+1), vn(nx+1,ny), pn(nx+1,ny+1)
        real(8) :: D(nx,ny)

        itc = itc+1
        error = 0.0d0

        do i=2,nx
            do j=2,ny
                D(i,j) = (ur(i,j)-ur(i-1,j))/dx+(vr(i,j)-vr(i,j-1))/dy
                if(D(i,j).GT.error) error = D(i,j)
            enddo
        enddo

        u = un
        v = vn
        pGuess = pn

        write(*,*) itc,' ',error

!!!        open(unit=01,file='error.dat',status='unknown',position='append')
!!!        if (MOD(itc,100).EQ.0) then
!!!            write(01,*) itc,' ',error
!!!        endif
!!!        close(01)

        return
        end subroutine check


!!! compute velocity components u, v and pressure p
        subroutine caluvp(uc,vc,pc)
        implicit none
        include "para.h"
        include "common.h"
        integer :: i, j
        real(8) :: uc(nx,ny), vc(nx,ny), pc(nx,ny)

        do i=1,nx
            do j=1,ny
                uc(i,j) = 0.5d0*(u(i,j)+u(i,j+1))
                vc(i,j) = 0.5d0*(v(i,j)+v(i+1,j))
                pc(i,j) = 0.25d0*(p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))
            enddo
        enddo

        return
        end subroutine caluvp


!!! compute Streamfunction
        subroutine calpsi(uc,vc)
        implicit none
        include "para.h"
        include "common.h"
        integer :: i, j
        real(8) :: uc(nx,ny), vc(nx,ny)

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
            psi(i,j+1) = uc(i,j)*2.0d0*dy+psi(i,j-1)
            !psi(i+1,j) = -v(i-1,j)*2.0d0*dx+psi(i-1,j) ! Alternative and equivalent psi formulae
            enddo
        enddo

        do j=2,ny-1
            psi(2,j) = 0.25d0*psi(3,j)
            psi(nx-1,j) = 0.25d0*psi(nx-2,j)
        enddo
        do i=2,nx-1
            psi(i,2) = 0.25d0*psi(i,3)
            psi(i,ny-1) = 0.25d0*(psi(i,ny-2)-2.0d0*dy)
        enddo

        return
        end subroutine calpsi

!!! output data file
    subroutine output(uc,vc)
    implicit none
    include "para.h"
    include "common.h"
    integer :: i, j
    real(8) :: uc(nx,ny), vc(nx,ny)
    character (len=100):: filename

    write(filename,*) itc
    filename = adjustl(filename)
    open(unit=02,file='lid_'//trim(filename)//'.plt',status='unknown')

    write(02,*) 'TITLE="Lid Driven Cavity Flow"'
    write(02,*) 'VARIABLES=x,y,u,v,psi'
    write(02,*) 'ZONE','i=',nx,'j=',ny,'F=POINT'

    do j=1,ny
        do i=1,nx
            write(02,100) X(i), Y(j), uc(i,j), vc(i,j), psi(i,j)
        enddo
    enddo

100 format(2x,10(e12.6,' '))

    close(02)

    return
    end subroutine output
