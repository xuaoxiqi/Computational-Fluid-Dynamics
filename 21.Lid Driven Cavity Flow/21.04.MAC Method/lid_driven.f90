
!!!    This program sloves Lid Driven Cavity Flow problem using Marker in Cell Method
!!!    Copyright (C) 2012  Ao Xu
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
        integer, parameter :: N=81,M=81
        integer :: itc, itc_max, k
        real(8) :: u(N,M+1),v(N+1,M),p(N+1,M+1),psi(N,M),X(N), Y(M)
        real(8) :: un(N,M+1),vn(N+1,M),pn(N+1,M+1),uc(N,M),vc(N,M),pc(N,M)
        real(8) :: Re, dt, dx, dy, eps, error

!!! input initial data
        Re = 1000.0d0
        dt = 1e-4
        dx = 1.0d0/float(N-1)
        dy = 1.0d0/float(M-1)
        eps = 1e-4
        itc = 0
        itc_max = 5*1e5
        error=100.00d0
        k = 0

!!! set up initial flow field
        call initial(N,M,dx,dy,X,Y,u,v,p,psi)

        do while((error.GT.eps).AND.(itc.LT.itc_max))

            error=0.0d0

!!! Solve Pressure Equation
            call calpn(N,M,dx,dy,dt,Re,p,u,v,pn)

!!! Solve Momentum Equation
            call solmom(N,M,dx,dy,dt,Re,u,v,pn,un,vn)

!!! check convergence
            call check(N,M,dx,dy,dt,error,u,v,p,un,vn,pn,itc)

!!! output preliminary results
            if (MOD(itc,10000).EQ.0) then

                !!! compute velocity components u, v and pressure p
                call caluvp(N,M,u,v,p,uc,vc,pc)

                call calpsi(N,M,dx,dy,uc,vc,psi)
                k = k+1
                call output(N,M,X,Y,uc,vc,psi,k)
            endif

        enddo
!!! compute velocity components u, v and pressure p
        call caluvp(N,M,u,v,p,uc,vc,pc)

!!! compute Streamfunction
        call calpsi(N,M,dx,dy,uc,vc,psi)

!!! output data file
        call output(N,M,X,Y,uc,vc,psi,k)

        write(*,*)
        write(*,*) '************************************************************'
        write(*,*) 'This program sloves Lid Driven Cavity Flow problem'
        write(*,*) 'using Marker in Cell Method'
        write(*,*) 'N =',N,',       M =',M
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
        subroutine initial(N,M,dx,dy,X,Y,u,v,p,psi)
        implicit none
        integer :: N, M, i, j
        real(8) :: dx, dy
        real(8) :: u(N,M+1), v(N+1,M), p(N+1,M+1), psi(N,M), uc(N,M), vc(N,M), X(N), Y(M)

        do i=1,N
            X(i) = (i-1)*dx
        enddo
        do j=1,M
            Y(j) = (j-1)*dy
        enddo

        p = 1.0d0
        u = 0.0d0
        v = 0.0d0
        psi = 0.0d0

        do i=1,N
            u(i,M+1) = 4.0d0/3.0d0
            u(i,M) = 2.0d0/3.0d0
        enddo

        return
        end subroutine initial

!!! Solve Pressure Equation
        subroutine calpn(N,M,dx,dy,dt,Re,p,u,v,pn)
        implicit none
        integer :: N, M, i, j
        real(8) :: Re
        real(8) :: dx, dy, dt
        real(8) :: alpha
        real(8) :: u(N,M+1), v(N+1,M), p(N+1,M+1), pn(N+1,M+1), D(N+1,M+1)
        real(8) :: R(N,M)

        D = 0.0d0
        do i=2,N
            do j=2,M
                D(i,j) = (u(i,j)-u(i-1,j))/dx+(v(i,j)-v(i,j-1))/dy
            enddo
        enddo

        do i=2,N
            do j=2,M
        R(i,j) = -(0.0d0-D(i,j))/dt -  (u(i,j)*u(i,j)-(u(i,j)+u(i-1,j))*(u(i,j)+u(i-1,j))/2.0d0+u(i-1,j)*u(i-1,j)) /dx/dx &
                - 2.0*(  0.25d0*( (u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j))+(u(i-1,j)+u(i-1,j-1))*(v(i-1,j-1)+v(i,j-1)) &
                                 -(u(i,j)+u(i,j-1))*(v(i,j-1)+v(i+1,j-1))+(u(i-1,j)+u(i-1,j+1))*(v(i-1,j)+v(i,j))  )  )/dx/dy &
                - (v(i,j)*v(i,j)-(v(i,j)+v(i,j-1))*(v(i,j)+v(i,j-1))/2.0d0+v(i,j-1)*v(i,j-1)) &
                + 1.0d0/Re*( (D(i+1,j)-2.0d0*D(i,j)+D(i-1,j))/dx/dx+(D(i,j+1)-2.0d0*D(i,j)+D(i,j-1))/dy/dy )
            enddo
        enddo

        alpha =0.6d0    !!! alpha is relaxation factor

        do i=2,N
            do j=2,M
                pn(i,j) = (1-alpha)*p(i,j)+alpha/2.0d0/(1.0d0/dx/dx+1.0d0/dy/dy)  &
                        *(1.0d0/dx/dx*(p(i+1,j)+p(i-1,j))+1.0d0/dy/dy*(p(i,j+1)+p(i,j-1))-R(i,j))
            enddo
        enddo

        !!! boundary condition
        do i=2,N
            pn(i,1) = pn(i,2)
            pn(i,M+1) = pn(i,M)
        enddo
        do j=1,M+1
            pn(1,j) = pn(2,j)
            pn(N+1,j) = pn(N,j)
        enddo

        return
        end subroutine calpn


!!! Solve Momentum Equation
        subroutine solmom(N,M,dx,dy,dt,Re,u,v,pn,un,vn)
        implicit none
        integer :: N, M, i, j
        real(8) :: u(N,M+1),v(N+1,M),pn(N+1,M+1),un(N,M+1),vn(N+1,M)
        real(8) :: Re, dx, dy, dt

        do i=2,N-1
            do j=2,M
    un(i,j) = u(i,j) - dt*(  (u(i+1,j)*u(i+1,j)-u(i-1,j)*u(i-1,j))/2.0d0/dx &
    +0.25d0*( (u(i,j)+u(i,j+1))*(v(i,j)+v(i-1,j))-(u(i,j)+u(i,j-1))*(v(i-1,j-1)+v(i,j-1)) )/dy  )&
    - dt/dx*(pn(i+1,j)-pn(i,j)) &
    + dt*1.0d0/Re*( (u(i+1,j)-2.0d0*u(i,j)+u(i-1,j))/dx/dx +(u(i,j+1)-2.0d0*u(i,j)+u(i,j-1))/dy/dy )
            enddo
        enddo

        do j=2,M
            un(1,j) = 0.0d0
            un(N,j) = 0.0d0
        enddo

        do i=1,N
            un(i,1) = -un(i,2)
            un(i,M+1) = 2.0d0-un(i,M)
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=2,N
            do j=2,M-1
    vn(i,j) = v(i,j) - dt* ( 0.25d0*( (u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j))-(u(i-1,j)+u(i-1,j+1))*(v(i,j)+v(i-1,j)) )/dx &
    +(v(i,j+1)*v(i,j+1)-v(i,j-1)*v(i,j-1))/2.0d0/dy ) &
    - dt/dy*(pn(i,j+1)-pn(i,j)) &
    + dt*1.0d0/Re*( (v(i+1,j)-2.0d0*v(i,j)+v(i-1,j))/dx/dx+(v(i,j+1)-2.0d0*v(i,j)+v(i,j-1))/dy/dy )
            enddo
        enddo

        do j=2,M-1
            vn(1,j) = -vn(2,j)
            vn(M+1,j) = -vn(M,j)
        enddo

        do i=1,N+1
            vn(i,1) = 0.0d0
            vn(i,M) = 0.0d0
        enddo

        return
        end subroutine solmom


!!! check convergence
        subroutine check(N,M,dx,dy,dt,error,u,v,p,un,vn,pn,itc)
        implicit none
        integer :: N, M, i, j, itc
        real(8) :: dt, error, temp
        real(8) :: dx, dy
        real(8) :: u(N,M+1), v(N+1,M), p(N+1,M+1), un(N,M+1), vn(N+1,M), pn(N+1,M+1), Dn(N,M)
        real(8) :: erru, errv, errp

        itc = itc+1
        error = 0.0d0
        Dn = 0.0d0

        do i=2,N
            do j=2,M
                Dn(i,j) = (un(i,j)-un(i-1,j))/dx+(vn(i,j)-vn(i,j-1))/dy
                if(Dn(i,j).GT.error) error = Dn(i,j)
            enddo
        enddo

        u = un
        v = vn
        p = pn


        write(*,*) itc,' ',error

!!!        open(unit=01,file='error.dat',status='unknown',position='append')
!!!        if (MOD(itc,100).EQ.0) then
!!!            write(01,*) itc,' ',error
!!!        endif
!!!        close(01)

        return
        end subroutine check


!!! compute velocity components u, v and pressure p
        subroutine caluvp(N,M,u,v,p,uc,vc,pc)
        implicit none
        integer :: N, M, i, j
        real(8) :: u(N,M+1), v(N+1,M), p(N+1,M+1), uc(N,M), vc(N,M), pc(N,M)

        do i=1,N
            do j=1,M
                uc(i,j) = 0.5d0*(u(i,j)+u(i,j+1))
                vc(i,j) = 0.5d0*(v(i,j)+v(i+1,j))
                pc(i,j) = 0.25d0*(p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))
            enddo
        enddo

        return
        end subroutine caluvp


!!! compute Streamfunction
        subroutine calpsi(N,M,dx,dy,u,v,psi)
        implicit none
        integer :: N, M, i, j
        real(8) :: dx, dy
        real(8) :: u(N,M), v(N,M), psi(N,M)

!        do j=1,M
!            psi(1,j) = 0.0d0
!            psi(N,j) = 0.0d0
!        enddo
!        do i=1,N
!            psi(i,1) = 0.0d0
!            psi(i,M) = 0.0d0
!        enddo

        do i=3,N-2
            do j=2,M-3
            psi(i,j+1) = u(i,j)*2.0d0*dy+psi(i,j-1)
            !psi(i+1,j) = -v(i-1,j)*2.0d0*dx+psi(i-1,j) ! Alternative and equivalent psi formulae
            enddo
        enddo

        do j=2,M-1
            psi(2,j) = 0.25d0*psi(3,j)
            psi(N-1,j) = 0.25d0*psi(N-2,j)
        enddo
        do i=2,N-1
            psi(i,2) = 0.25d0*psi(i,3)
            psi(i,M-1) = 0.25d0*(psi(i,M-2)-2.0d0*dy)
        enddo


        return
        end subroutine calpsi

!!! output data file
        subroutine output(N,M,X,Y,uc,vc,psi,k)
        implicit none
        integer :: N, M, i, j, k
        real(8) :: X(N), Y(M), uc(N,M), vc(N,M), psi(N,M)

        character*16 filename

        filename='0000cavity.dat'
        filename(1:1) = CHAR(ICHAR('0')+MOD(k/1000,10))
        filename(2:2) = CHAR(ICHAR('0')+MOD(k/100,10))
        filename(3:3) = CHAR(ICHAR('0')+MOD(k/10,10))
        filename(4:4) = CHAR(ICHAR('0')+MOD(k,10))

        open(unit=02,file=filename,status='unknown')
        write(02,101)
        write(02,102)
        write(02,103) N, M
        do j=1,M
            do i = 1,N
                write(02,100) X(i), Y(j), uc(i,j), vc(i,j), psi(i,j)
            enddo
        enddo

100     format(2x,10(e12.6,'      '))
101     format('Title="Lid Driven Cavity Flow"')
102     format('Variables=x,y,u,v,psi')
103     format('zone',1x,'i=',1x,i5,2x,'j=',1x,i5,1x,'f=point')

        close(02)

        return
        end subroutine output



