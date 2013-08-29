

!!!    This program sloves Buoyancy Driven Cavity Flow problem using Vorticity-Streamfunction Methods
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

!!!           Adiabatic Wall
!!!         |---------------|
!!!         |               |
!!!         |       |       |
!!!    Hot  |       |       | Cold
!!!    Wall |       | g     | Wall
!!!         |       v       |
!!!         |               |
!!!         |---------------|
!!!           Adiabatic Wall

        program main
        implicit none
        integer, parameter :: nx=81, ny=81
        integer :: i, j, itc, itc_max, k
        real(8) :: dx, dy, Pr, Ra, dt, eps, error
        real(8) :: X(nx), Y(ny), u(nx,ny), v(nx,ny), vor(nx,ny), RVOR(nx,ny), psi(nx,ny), Rpsi(nx,ny), T(nx,ny)

!!! input initial data
        Pr = 0.7d0
        Ra = 1e3
        dx = 1.0d0/(nx-1)
        dy = 1.0d0/(ny-1)
        dt = 3*1e-5
        eps = 1e-7
        itc = 0
        itc_max = 1e8
        error = 100.0d0

!!! set up initial flow field
        call initial(nx,ny,dx,dy,X,Y,u,v,psi,vor,RVOR,Rpsi,T)

        do while((error.GT.eps).AND.(itc.LT.itc_max))
!!! solve vorticity equation
            call solvor(nx,ny,dx,dy,Pr,Ra,dt,u,v,vor,RVOR,T,psi)

!!! solve Streamfunction equation
            call solpsi(nx,ny,dx,dy,vor,psi,Rpsi)

!!! compute velocity components u and v
            call caluv(nx,ny,dx,dy,psi,u,v)

!!! compute temperature field
            call calT(nx,ny,dt,dx,dy,u,v,T)

!!! check convergence
            call convergence(nx,ny,dt,RVOR,Rpsi,error,itc)

!!! output preliminary results
            if (MOD(itc,10000).EQ.0) then
                call output(nx,ny,X,Y,u,v,psi,T,itc)
            endif

        enddo

!!! output data file
        call output(nx,ny,X,Y,u,v,psi,T,itc)

        open(unit=03,file='results.log',status='unknown')

        write(03,*)
        write(03,*) '************************************************************'
        write(03,*) 'This program sloves Buoyancy Driven Cavity Flow problem'
        write(03,*) 'using Vorticity-Streamfunction Methods'
        write(03,*) 'nx =',nx,',       ny =',ny
        write(03,*) 'Pr=',Pr
        write(03,*) 'Ra=',Ra
        write(03,*) 'dt =', dt
        write(03,*) 'eps =',eps
        write(03,*) 'itc =',itc
        write(03,*) 'Developing time=',dt*itc,'s'
        write(03,*) '************************************************************'
        write(03,*)

        close(03)

        stop
        end program main


!!! set up initial flow field
        subroutine initial(nx,ny,dx,dy,X,Y,u,v,psi,vor,RVOR,Rpsi,T)
        implicit none
        integer :: i, j, nx, ny
        real(8) :: dx, dy
        real(8) :: X(nx), Y(ny), u(nx,ny), v(nx,ny), psi(nx,ny), vor(nx,ny), RVOR(nx,ny), Rpsi(nx,ny), T(nx,ny)

        do i=1,nx
            X(i) = (i-1)*dx
        enddo
        do j=1,ny
            Y(j) = (j-1)*dy
        enddo

        u = 0.0d0
        v = 0.0d0
        psi = 0.0d0
        Rpsi = 0.00d0
        vor = 0.0d0
        RVOR = 0.0d0
        T = 0.0d0

        do j =1,ny
            T(1,j) = 1.0d0  !Left side hot wall
        enddo

        return
        end subroutine initial


!!! solve vorticity equation
        subroutine solvor(nx,ny,dx,dy,Pr,Ra,dt,u,v,vor,RVOR,T,psi)
        implicit none
        integer :: i, j, nx, ny
        real(8) :: dx, dy, dt, Pr, Ra, dvorx2, dvory2, dvorx1, dvory1
        real(8) :: vor(nx,ny), u(nx,ny), v(nx,ny), RVOR(nx,ny), T(nx,ny), psi(nx,ny)
        real(8) :: vori(nx,ny)

        call RESvor(nx,ny,dx,dy,Pr,Ra,dt,u,v,vor,RVOR,T)
        do i=2,nx-1
            do j=2,ny-1
                vori(i,j) = vor(i,j)+dt*RVOR(i,j)
            enddo
        enddo
        call bcvor(nx,ny,dx,dy,vori,psi)

        call RESvor(nx,ny,dx,dy,Pr,Ra,dt,u,v,vori,RVOR,T)
        do i=2,nx-1
            do j=2,ny-1
                vori(i,j) = 0.75d0*vor(i,j)+0.25d0*(vori(i,j)+dt*RVOR(i,j))
            enddo
        enddo
        call bcvor(nx,ny,dx,dy,vori,psi)

        call RESvor(nx,ny,dx,dy,Pr,Ra,dt,u,v,vori,RVOR,T)
        do i=2,nx-1
            do j=2,ny-1
                vor(i,j) = 1.0d0/3.0d0*vor(i,j)+2.0d0/3.0d0*(vori(i,j)+dt*RVOR(i,j))
            enddo
        enddo
        call bcvor(nx,ny,dx,dy,vor,psi)

        return
        end subroutine solvor

        subroutine RESvor(nx,ny,dx,dy,Pr,Ra,dt,u,v,vor,RVOR,T)
        implicit none
        integer :: i, j, nx, ny
        real(8) :: dx, dy, dt, Pr, Ra, dvorx2, dvory2, dvorx1, dvory1
        real(8) :: vor(nx,ny), u(nx,ny), v(nx,ny), RVOR(nx,ny), T(nx,ny)

        do i=2,nx-1
            do j=2,ny-1
                dvorx2 = (vor(i+1,j)-2.0d0*vor(i,j)+vor(i-1,j))/dx/dx
                dvory2 = (vor(i,j+1)-2.0d0*vor(i,j)+vor(i,j-1))/dy/dy
                dvorx1 = (u(i+1,j)*vor(i+1,j)-u(i-1,j)*vor(i-1,j))/2.0d0/dx
                dvory1 = (v(i,j+1)*vor(i,j+1)-v(i,j-1)*vor(i,j-1))/2.0d0/dy
                RVOR(i,j) = (dvorx2+dvory2)*Pr-Ra*Pr*(T(i+1,j)-T(i-1,j))/2.0d0/dx-dvorx1-dvory1
            enddo
        enddo

        return
        end subroutine RESvor


!!! updates the boundary condition for vorticity
        subroutine bcvor(nx,ny,dx,dy,vor,psi)
        implicit none
        integer :: i, j, nx, ny
        real(8) :: dx, dy
        real(8) :: vor(nx,ny), psi(nx,ny)

        ! 2nd order approximation
        do j=1,ny
            vor(1,j) = 3.0d0*psi(2,j)/dx/dx-0.5d0*vor(2,j)
            vor(nx,j) = 3.0d0*psi(nx-1,j)/dx/dx-0.5d0*vor(nx-1,j)
        enddo
        do i=1,nx
            vor(i,1) = 3.0d0*psi(i,2)/dy/dy-0.5d0*vor(i,2)
            vor(i,ny) = 3.0d0*psi(i,ny-1)/dy/dy-0.5d0*vor(i,ny-1)
        enddo

        return
        end subroutine bcvor

!!! solve Streamfunction equation
        subroutine solpsi(nx,ny,dx,dy,vor,psi,Rpsi)
        implicit none
        integer :: i, j ,nx, ny
        real(8) :: alpha, dx, dy, aw, as, ap
        real(8) :: vor(nx,ny), psi(nx,ny), Rpsi(nx,ny), S(nx,ny)

        aw = 1.0d0/dx/dx
        as = 1.0d0/dy/dy
        ap = -2.0d0*(as+aw)

        do i=3,nx-2
            do j=3,ny-2
                S(i,j) = vor(i,j)-(psi(i+1,j)-2.0d0*psi(i,j)+psi(i-1,j))/dx/dx &
                                -(psi(i,j+1)-2.0d0*psi(i,j)+psi(i,j-1))/dy/dy
            enddo
        enddo

        do j=1,ny
            Rpsi(1,j) = 0.0d0
            Rpsi(2,j) = 0.0d0
            Rpsi(nx,j) = 0.0d0
            Rpsi(nx-1,j) = 0.0d0
        enddo
        do i=1,nx
            Rpsi(i,1) = 0.0d0
            Rpsi(i,2) = 0.0d0
            Rpsi(i,ny) = 0.0d0
            Rpsi(i,ny-1) = 0.0d0
        enddo

        alpha = 1.3d0      !alpha is ralaxtion factor

        do i=3,nx-2
            do j=3,ny-2
                Rpsi(i,j)=(S(i,j)-aw*Rpsi(i-1,j)-as*Rpsi(i,j-1))/ap
                psi(i,j) = psi(i,j)+alpha*Rpsi(i,j)
            enddo
        enddo


        call bcpsi(nx,ny,dy,psi)

        return
        end subroutine solpsi


!!! updates the values of sream function at boundary points
        subroutine bcpsi(nx,ny,dy,psi)
        implicit none
        integer :: i, j, nx, ny
        real(8) :: dy
        real(8) :: psi(nx,ny)

        do j=2,ny-1
            psi(2,j) = 0.25d0*psi(3,j)
            psi(nx-1,j) = 0.25d0*psi(nx-2,j)
        enddo
        do i=2,nx-1
            psi(i,2) = 0.25d0*psi(i,3)
            psi(i,ny-1) = 0.25d0*psi(i,ny-2)
        enddo


        return
        end subroutine bcpsi

!!! compute velocity components u and v
        subroutine caluv(nx,ny,dx,dy,psi,u,v)
        implicit none
        integer :: i, j, nx, ny
        real(8) :: dx, dy
        real(8) :: psi(nx,ny), u(nx,ny), v(nx,ny)
        integer velocity_order

        !physical boundary condition
        do i=1,nx
            u(i,1) = 0.0d0
            v(i,1) = 0.0d0
            u(i,ny) = 0.0d0
            v(i,ny) = 0.0d0
        enddo
        do j=1,ny
            u(1,j) = 0.0d0
            u(nx,j) = 0.0d0
            v(1,j) = 0.0d0
            v(nx,j) = 0.0d0
        enddo


        do i=2,nx-1
            do j=2,ny-1
                u(i,j) = 0.5d0*(psi(i,j+1)-psi(i,j-1))/dy
                v(i,j) = -0.5d0*(psi(i+1,j)-psi(i-1,j))/dx
            enddo
        enddo

        return
        end subroutine caluv

!!! compute temperature field
        subroutine calT(nx,ny,dt,dx,dy,u,v,T)
        implicit none
        integer :: i, j, nx, ny
        real(8) :: dx, dy, dt, dTx2, dTy2, dTx1, dTy1
        real(8) :: T(nx,ny), u(nx,ny), v(nx,ny), RT(nx,ny), Ti(nx,ny)

        call REST(nx,ny,dx,dy,dt,u,v,T,RT)
        do i=2,nx-1
            do j=2,ny-1
                Ti(i,j) = T(i,j)+dt*RT(i,j)
            enddo
        enddo
        call bcT(nx,ny,dx,dy,dt,Ti)

        call REST(nx,ny,dx,dy,dt,u,v,Ti,RT)
        do i=2,nx-1
            do j=2,ny-1
                Ti(i,j) = 0.75d0*T(i,j)+0.25d0*(Ti(i,j)+dt*RT(i,j))
            enddo
        enddo
        call bcT(nx,ny,dx,dy,dt,Ti)

        call REST(nx,ny,dx,dy,dt,u,v,Ti,RT)
        do i=2,nx-1
            do j=2,ny-1
                T(i,j) = 1.0d0/3.0d0*T(i,j)+2.0d0/3.0d0*(Ti(i,j)+dt*RT(i,j))
            enddo
        enddo
        call bcT(nx,ny,dx,dy,dt,T)

        return
        end subroutine calT


        subroutine REST(nx,ny,dx,dy,dt,u,v,T,RT)
        implicit none
        integer :: i, j, nx, ny
        real(8) :: dx, dy, dt, dTx2, dTy2, dTx1, dTy1
        real(8) :: T(nx,ny), u(nx,ny), v(nx,ny), RT(nx,ny)


        ! Interior points using FTCS Sheme
        do i=2,nx-1
            do j=2,ny-1
                dTx2 = (T(i+1,j)-2.0d0*T(i,j)+T(i-1,j))/dx/dx
                dTy2 = (T(i,j+1)-2.0d0*T(i,j)+T(i,j-1))/dy/dy
                dTx1 = (u(i+1,j)*T(i+1,j)-u(i-1,j)*T(i-1,j))/2.0/dx
                dTy1 = (v(i,j+1)*T(i,j+1)-v(i,j-1)*T(i,j-1))/2.0/dy
                RT(i,j) = dTx2+dTy2-dTx1-dTy1
            enddo
        enddo

        return
        end subroutine REST


        subroutine bcT(nx,ny,dx,dy,dt,T)
        implicit none
        integer :: i, j, nx, ny
        real(8) :: dx, dy, dt
        real(8) :: T(nx,ny)
        integer :: bc_order

        !Left and right side boundary(Dirichlet B.C.)
        do j=1,ny
            T(1,j) = 1.0d0
            T(nx,j) = 0.0d0
        enddo

        do i=1,nx
            T(i,1) = (4.0d0*T(i,2)-T(i,3))/3.0d0
            T(i,ny) = (4.0d0*T(i,ny-1)-T(i,ny-2))/3.0d0
        enddo

        return
        end subroutine bcT


!!! check convergence
        subroutine convergence(nx,ny,dt,RVOR,Rpsi,error,itc)
        implicit none
        integer :: nx, ny, i, j, itc
        real(8) :: RVOR(nx,ny), Rpsi(nx,ny)
        real(8) :: dt, error, errvor, errpsi

        itc = itc+1
        errvor = 0.0d0
        errpsi = 0.0d0

        do i=1,nx
            do j=1,ny
                if(ABS(RVOR(i,j))*dt.GT.errvor) errvor = ABS(RVOR(i,j))*dt
                if(ABS(Rpsi(i,j)).GT.errpsi) errpsi = ABS(Rpsi(i,j))
            enddo
        enddo

        error = MAX(errvor,errpsi)
        if(itc.EQ.1) error = 100.0d0



        open(unit=01,file='error.dat',status='unknown',position='append')

        ! write(*,*) itc,' ',error
        if (MOD(itc,2000).EQ.0) then
            write(*,*) itc,' ',error
            write(01,*) itc,' ',error
        endif

        close(01)


        return
        end subroutine convergence


!!! output data file
        subroutine output(nx,ny,X,Y,u,v,psi,T,itc)
        implicit none
        integer :: nx, ny, i, j, itc
        real(8) :: X(nx), Y(ny), u(nx,ny), v(nx,ny), psi(nx,ny),T(nx,ny)
        character(len=100) :: filename

        write(filename,*) itc
        filename = adjustl(filename)

        open(unit=02,file='cavity_'//trim(filename)//'.plt',status='unknown')
        write(02,*) 'TITLE="Buoyancy Driven Cavity Flow"'
        write(02,*) 'VARIABLES=x,y,u,v,psi,T'
        write(02,101) nx, ny
        do j=1,ny
            do i = 1,nx
                write(02,100) X(i), Y(j), u(i,j), v(i,j), psi(i,j), T(i,j)
            enddo
        enddo

100     format(2x,10(e12.6,' '))
101     format('ZONE',1x,'i=',1x,i5,2x,'j=',1x,i5,1x,'F=POINT')

        close(02)

        return
        end subroutine output



