
!!!    This program sloves Lid Driven Cavity Flow problem using Lattice Boltzmann Method
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
        integer :: i, j, itc, itc_max, k
        integer :: iwall(N,M)
        real(8) :: Re, cs2, U_ref, L_ref, dx, dy, tau, omega
        real(8) :: eps, error
        real(8) :: X(N), Y(M), u(N,M), v(N,M), up(N,M), vp(N,M), rho(N,M), p(N,M), psi(N,M)
        real(8) :: t_k(0:8), f(0:8,N,M), un(0:8)
        real(8) :: xv(9), yv(9)
        data xv/0.0d0,1.0d0,0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0/
        data yv/0.0d0,0.0d0,1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0/

!!! input initial data
        Re = 100.0d0
        cs2 = 1.0d0/3.0d0
        U_ref = 1.0d0
        L_ref = 1.0d0
        dx = L_ref/N
        dy = L_ref/M
        do i=1,N
            X(i) = (i-1)*dx
        enddo
        do j=1,M
            Y(j) = (j-1)*dy
        enddo
        tau =  1.0d0/cs2*U_ref*(L_ref-1.d0)/Re+0.5d0
        omega = 1.0d0/tau
        itc = 0
        itc_max = 1000
        eps = 1e-8
        k = 0
!!! set up initial flow field
        call initial(N,M,u,v,rho,psi,iwall,U_ref,cs2,t_k,xv,yv,un,f)

        do while(itc.LT.itc_max)

            call propagate(N,M,f)

            call relaxation(N,M,iwall,u,v,xv,yv,rho,f,t_k,cs2,omega)

            call bounceback(N,M,f)

            call caluv(N,M,u,v,up,vp)

            call check(N,M,iwall,u,v,up,vp,itc,error)

!!! output preliminary results
            if (MOD(itc,10000).EQ.0) then

!                !!! compute velocity components u, v and pressure p
!                call caluvp(N,M,u,v,p,uc,vc,pc)
!
!                call calpsi(N,M,dx,dy,uc,vc,psi)
!                k = k+1
!                call output(N,M,X,Y,uc,vc,psi,k)
            endif

        enddo

!!! compute Streamfunction
        call calpsi(N,M,dx,dy,up,vp,psi)

!!! output data file
        call output(N,M,X,Y,up,vp,psi,k)

        write(*,*)
        write(*,*) '************************************************************'
        write(*,*) 'This program sloves Lid Driven Cavity Flow problem'
        write(*,*) 'using Lattice Boltzmann Method'
        write(*,*) 'Consider D2Q9 Model'
        write(*,*) 'N =',N,',       M =',M
        write(*,*) 'Re =',Re
        write(*,*) 'eps =',eps
        write(*,*) 'itc =',itc
        write(*,*) '************************************************************'
        write(*,*)

        stop
        end program main

        subroutine initial(N,M,u,v,rho,psi,iwall,U_ref,cs2,t_k,xv,yv,un,f)
        implicit none
        integer :: x, y, N, M, i
        integer :: iwall(N,M)
        real(8) :: U_ref, cs2, us2
        real(8) :: t_k(0:8), u(N,M), v(N,M), rho(N,M), psi(N,M), xv(0:8), yv(0:8), un(0:8)
        real(8) :: f(0:8,N,M)

        psi = 0.0d0

        t_k(0) = 4.0d0/9.0d0
        do i=1,4
            t_k(i) = 1.0d0/9.0d0
        enddo
        do i=5,8
            t_k(i) = 1.0d0/36.0d0
        enddo

        do x=1,N
            do y=1,M
                iwall(x,y) = 0
                u(x,y) = 0.0d0
                if(y.EQ.M) u(x,y) = U_ref
                v(x,y) = 0.0d0
                rho(x,y) = 1.0d0
            enddo
        enddo

        !Wall type
        do x=1,N
            iwall(x,1) = 1
            iwall(x,M) = 2
        enddo
        do y=1,M
            iwall(1,y) = 1
            iwall(N,y) = 1
        enddo

        do x=1,N
            do y=1,M
                us2 = u(x,y)*u(x,y)+v(x,y)*v(x,y)
                do i=0,8
                    un(i) = u(x,y)*xv(i)+v(x,y)*yv(i)
                    f(i,x,y) = t_k(i)*(1.0d0+un(i)/cs2+un(i)*un(i)/(2.0d0*cs2*cs2)-us2/(2.0d0*cs2))
                enddo
            enddo
        enddo

        return
        end subroutine initial

        subroutine propagate(N,M,f)
        implicit none
        integer :: x, y, N, M
        real(8) :: f(0:8,N,M)

        do x=1,N
            do y=1,M-1
                f(0,x,y) = f(0,x,y)
            enddo
        enddo
        do x=N,2,-1
            do y=1,M-1
                f(1,x,y) = f(1,x-1,y)
            enddo
        enddo
        do x=1,N
            do y=M-1,2,-1
                f(2,x,y) = f(2,x,y-1)
            enddo
        enddo
        do x=1,N-1
            do y=1,M-1
                f(3,x,y) = f(3,x+1,y)
            enddo
        enddo
        do x=1,N
            do y=1,M-1
                f(4,x,y) = f(4,x,y+1)
            enddo
        enddo
        do x=N,2,-1
            do y=M-1,2,-1
                f(5,x,y) = f(5,x-1,y-1)
            enddo
        enddo
        do x=1,N-1
            do y=M-1,2,-1
                f(6,x,y) = f(6,x+1,y-1)
            enddo
        enddo
        do x=1,N-1
            do y=1,M-1
                f(7,x,y) = f(7,x+1,y+1)
            enddo
        enddo
        do x=N,2,-1
            do y=1,M-1
                f(8,x,y) = f(8,x-1,y+1)
            enddo
        enddo

        return
        end subroutine propagate

        subroutine relaxation(N,M,iwall,u,v,xv,yv,rho,f,t_k,cs2,omega)
        implicit none
        integer :: x, y, N, M, i
        integer :: iwall(N,M)
        real(8) :: cs2, us2, omega
        real(8) :: rho(N,M), u(N,M), v(N,M), u2(N,M)
        real(8) :: xv(0:8), yv(0:8), un(0:8), t_k(0:8), feq(0:8,N,M), f(0:8,N,M)

        do x=1,N
            do y=1,M-1
                if(iwall(x,y).NE.2) then
                    rho(x,y) = 0.0d0
                    do i=0,8
                        rho(x,y) = rho(x,y)+f(i,x,y)
                    enddo
                    !data xv/0.0d0,1.0d0,0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0/
                    !data yv/0.0d0,0.0d0,1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0/
                    u(x,y) = (f(1,x,y)-f(3,x,y)+f(5,x,y)-f(6,x,y)-f(7,x,y)+f(8,x,y))/rho(x,y)
                    v(x,y) = (f(2,x,y)-f(4,x,y)+f(5,x,y)+f(6,x,y)-f(7,x,y)-f(8,x,y))/rho(x,y)
                    u2(x,y) = u(x,y)*u(x,y)+v(x,y)
                    do i=0,8
                        un(i) = u(x,y)*xv(i) + v(x,y)*yv(i)
                        feq(i,x,y) = t_k(i)*rho(x,y)*(1.0d0+un(i)/cs2+un(i)*un(i)/(2.0d0*cs2*cs2)-us2/(2.0d0*cs2))
                        f(i,x,y) = f(i,x,y)-omega*(f(i,x,y)-feq(i,x,y))
                    enddo

                    if((x.EQ.1).AND.(y.EQ.1)) then
                        f(6,x,y) = feq(6,x,y)
                        f(8,x,y) = feq(8,x,y)
                    endif

                    if((x.EQ.N).AND.(y.EQ.1)) then
                        f(5,x,y) = feq(5,x,y)
                        f(7,x,y) = feq(7,x,y)
                    endif

                endif

            enddo
        enddo

        return
        end subroutine relaxation

        subroutine bounceback(N,M,f)
        implicit none
        integer :: x, y, N, M
        real(8) :: f(0:8,N,M)

        do x=1,N
            do y=1,M-1

                !Left side
                if((x.EQ.1).AND.(y.NE.1)) then
                    f(1,x,y) = f(3,x,y)
                    f(5,x,y) = f(7,x,y)
                    f(8,x,y) = f(6,x,y)
                endif

                !Right side
                if((x.EQ.N).AND.(y.NE.1)) then
                    f(3,x,y) = f(1,x,y)
                    f(6,x,y) = f(8,x,y)
                    f(7,x,y) = f(5,x,y)
                endif

                if(y.EQ.1) then
                    !Left-Bottom corner
                    if(x.EQ.1) then
                        f(1,x,y) = f(3,x,y)
                        f(2,x,y) = f(4,x,y)
                        f(5,x,y) = f(7,x,y)
                    !Right-Bottom corner
                    elseif(x.EQ.N) then
                        f(3,x,y) = f(1,x,y)
                        f(6,x,y) = f(8,x,y)
                        f(2,x,y) = f(4,x,y)
                    !Bottom side
                    else
                        f(2,x,y) = f(4,x,y)
                        f(5,x,y) = f(7,x,y)
                        f(6,x,y) = f(8,x,y)
                    endif
                endif

            enddo
        enddo

        return
        end subroutine bounceback

        subroutine caluv(N,M,u,v,up,vp)
        implicit none
        integer :: x, y, N, M
        real(8) :: u(N,M), v(N,M), up(N,M), vp(N,M)

        do x=1,N
            do y=1,M
                up(x,y) = u(x,y)
                vp(x,y) = v(x,y)
            enddo
        enddo

        return
        end subroutine caluv


        subroutine check(N,M,iwall,u,v,up,vp,itc,error)
        implicit none
        integer :: x, y, N, M, itc
        integer :: iwall(N,M)
        real(8) :: error
        real(8) :: u(N,M), v(N,M), up(N,M), vp(N,M)

        itc = itc+1
        error = 0.0d0

        do x=1,N
            do y=1,M-1
                if(iwall(x,y).NE.2) then
                    error  = error+SQRT((u(x,y)-up(x,y))*(u(x,y)-up(x,y))+(v(x,y)-vp(x,y))*(v(x,y)-vp(x,y))) &
                                    /SQRT((u(x,y)+0.01)*(u(x,y)+0.01)+v(x,y)*v(x,y))
                endif
            enddo
        enddo

        write(*,*) itc,' ',error

!!!        open(unit=01,file='error.dat',status='unknown',position='append')
!!!        if (MOD(itc,2000).EQ.0) then
!!!            write(01,*) itc,' ',error
!!!        endif
!!!        close(01)

        return
        end subroutine check


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
        subroutine output(N,M,X,Y,u,v,psi,k)
        implicit none
        integer :: N, M, i, j, k
        real(8) :: X(N), Y(M), u(N,M), v(N,M), psi(N,M)

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
            do i=1,N
                write(02,100) X(i), Y(j), u(i,j), v(i,j), psi(i,j)
            enddo
        enddo

100     format(2x,10(e12.6,'      '))
101     format('Title="Lid Driven Cavity Flow"')
102     format('Variables=x,y,u,v,psi')
103     format('zone',1x,'i=',1x,i5,2x,'j=',1x,i5,1x,'f=point')

        close(02)

        return
        end subroutine output
