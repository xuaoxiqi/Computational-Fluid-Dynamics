

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
        integer, parameter :: N=81, M=81
        integer :: i, j, itc, itc_max, k
        real(8) :: dx, dy, Pr, Ra, dt, eps, error
        real(8) :: X(N), Y(M), u(N,M), v(N,M), vor(N,M), RVOR(N,M), psi(N,M), Rpsi(N,M), T(N,M)

!!! input initial data
        Pr = 0.71d0
        Ra = 1e4
        dx = 1.0d0/(N-1)
        dy = 1.0d0/(M-1)
        dt = 3*1e-5
        eps = 1e-8
        itc = 0
        itc_max = 1e8
        error = 100.0d0
        k = 0

!!! set up initial flow field
        call initial(N,M,dx,dy,X,Y,u,v,psi,vor,RVOR,Rpsi,T)

        do while((error.GT.eps).AND.(itc.LT.itc_max))
!!! solve vorticity equation
            call solvor(N,M,dx,dy,Pr,Ra,dt,u,v,vor,RVOR,T,psi)

!!! solve Streamfunction equation
            call solpsi(N,M,dx,dy,vor,psi,Rpsi)

!!! compute velocity components u and v
            call caluv(N,M,dx,dy,psi,u,v)

!!! compute temperature field
            call calT(N,M,dt,dx,dy,u,v,T)

!!! check convergence
            call convergence(N,M,dt,RVOR,Rpsi,error,itc)

!!! output preliminary results
            if (MOD(itc,10000).EQ.0) then
                k = k+1
                call output(N,M,X,Y,u,v,psi,T,k)
            endif

        enddo

!!! output data file
        k = k+1
        call output(N,M,X,Y,u,v,psi,T,k)

        open(unit=03,file='results.txt',status='unknown')

        write(03,*)
        write(03,*) '************************************************************'
        write(03,*) 'This program sloves Buoyancy Driven Cavity Flow problem'
        write(03,*) 'using Vorticity-Streamfunction Methods'
        write(03,*) 'N =',N,',       M =',M
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
        subroutine initial(N,M,dx,dy,X,Y,u,v,psi,vor,RVOR,Rpsi,T)
        implicit none
        integer :: i, j, N, M
        real(8) :: dx, dy
        real(8) :: X(N), Y(M), u(N,M), v(N,M), psi(N,M), vor(N,M), RVOR(N,M), Rpsi(N,M), T(N,M)

        do i=1,N
            X(i) = (i-1)*dx
        enddo
        do j=1,M
            Y(j) = (j-1)*dy
        enddo

        u = 0.0d0
        v = 0.0d0
        psi = 0.0d0
        Rpsi = 0.00d0
        vor = 0.0d0
        RVOR = 0.0d0
        T = 0.0d0

        do j =1,M
            T(1,j) = 1.0d0  !Left side hot wall
        enddo

        return
        end subroutine initial


!!! solve vorticity equation
        subroutine solvor(N,M,dx,dy,Pr,Ra,dt,u,v,vor,RVOR,T,psi)
        implicit none
        integer :: i, j, N, M
        real(8) :: dx, dy, dt, Pr, Ra, dvorx2, dvory2, dvorx1, dvory1
        real(8) :: vor(N,M), u(N,M), v(N,M), RVOR(N,M), T(N,M), psi(N,M)
        real(8) :: vori(N,M)

        call RESvor(N,M,dx,dy,Pr,Ra,dt,u,v,vor,RVOR,T)
        do i=2,N-1
            do j=2,M-1
                vori(i,j) = vor(i,j)+dt*RVOR(i,j)
            enddo
        enddo
        call bcvor(N,M,dx,dy,vori,psi)

        call RESvor(N,M,dx,dy,Pr,Ra,dt,u,v,vori,RVOR,T)
        do i=2,N-1
            do j=2,M-1
                vori(i,j) = 0.75d0*vor(i,j)+0.25d0*(vori(i,j)+dt*RVOR(i,j))
            enddo
        enddo
        call bcvor(N,M,dx,dy,vori,psi)

        call RESvor(N,M,dx,dy,Pr,Ra,dt,u,v,vori,RVOR,T)
        do i=2,N-1
            do j=2,M-1
                vor(i,j) = 1.0d0/3.0d0*vor(i,j)+2.0d0/3.0d0*(vori(i,j)+dt*RVOR(i,j))
            enddo
        enddo
        call bcvor(N,M,dx,dy,vor,psi)

        return
        end subroutine solvor

        subroutine RESvor(N,M,dx,dy,Pr,Ra,dt,u,v,vor,RVOR,T)
        implicit none
        integer :: i, j, N, M
        real(8) :: dx, dy, dt, Pr, Ra, dvorx2, dvory2, dvorx1, dvory1
        real(8) :: vor(N,M), u(N,M), v(N,M), RVOR(N,M), T(N,M)

        do i=2,N-1
            do j=2,M-1
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
        subroutine bcvor(N,M,dx,dy,vor,psi)
        implicit none
        integer :: i, j, N, M
        real(8) :: dx, dy
        real(8) :: vor(N,M), psi(N,M)

        ! 2nd order approximation
        do j=1,M
            vor(1,j) = 3.0d0*psi(2,j)/dx/dx-0.5d0*vor(2,j)
            vor(N,j) = 3.0d0*psi(N-1,j)/dx/dx-0.5d0*vor(N-1,j)
        enddo
        do i=1,N
            vor(i,1) = 3.0d0*psi(i,2)/dy/dy-0.5d0*vor(i,2)
            vor(i,M) = 3.0d0*psi(i,M-1)/dy/dy-0.5d0*vor(i,M-1)
        enddo

        return
        end subroutine bcvor

!!! solve Streamfunction equation
        subroutine solpsi(N,M,dx,dy,vor,psi,Rpsi)
        implicit none
        integer :: i, j ,N, M
        real(8) :: alpha, dx, dy, aw, as, ap
        real(8) :: vor(N,M), psi(N,M), Rpsi(N,M), S(N,M)

        aw = 1.0d0/dx/dx
        as = 1.0d0/dy/dy
        ap = -2.0d0*(as+aw)

        do i=3,N-2
            do j=3,M-2
                S(i,j) = vor(i,j)-(psi(i+1,j)-2.0d0*psi(i,j)+psi(i-1,j))/dx/dx-(psi(i,j+1)-2.0d0*psi(i,j)+psi(i,j-1))/dy/dy
            enddo
        enddo

        do j=1,M
            Rpsi(1,j) = 0.0d0
            Rpsi(2,j) = 0.0d0
            Rpsi(N,j) = 0.0d0
            Rpsi(N-1,j) = 0.0d0
        enddo
        do i=1,N
            Rpsi(i,1) = 0.0d0
            Rpsi(i,2) = 0.0d0
            Rpsi(i,M) = 0.0d0
            Rpsi(i,M-1) = 0.0d0
        enddo

        alpha = 1.5d0      !alpha is ralaxtion factor

        do i=3,N-2
            do j=3,M-2
                Rpsi(i,j)=(S(i,j)-aw*Rpsi(i-1,j)-as*Rpsi(i,j-1))/ap
                psi(i,j) = psi(i,j)+alpha*Rpsi(i,j)
            enddo
        enddo

        call bcpsi(N,M,dy,psi)

        return
        end subroutine solpsi


!!! updates the values of sream function at boundary points
        subroutine bcpsi(N,M,dy,psi)
        implicit none
        integer :: i, j, N, M
        real(8) :: dy
        real(8) :: psi(N,M)

        do j=2,M-1
            psi(2,j) = 0.25d0*psi(3,j)
            psi(N-1,j) = 0.25d0*psi(N-2,j)
        enddo
        do i=2,N-1
            psi(i,2) = 0.25d0*psi(i,3)
            psi(i,M-1) = 0.25d0*psi(i,M-2)
        enddo


        return
        end subroutine bcpsi


!!! compute velocity components u and v
        subroutine caluv(N,M,dx,dy,psi,u,v)
        implicit none
        integer :: i, j, N, M
        real(8) :: dx, dy
        real(8) :: psi(N,M), u(N,M), v(N,M)

        !physical boundary condition
        do i=1,N
            u(i,1) = 0.0d0
            v(i,1) = 0.0d0
            u(i,M) = 0.0d0
            v(i,M) = 0.0d0
        enddo
        do j=1,M
            u(1,j) = 0.0d0
            u(N,j) = 0.0d0
            v(1,j) = 0.0d0
            v(N,j) = 0.0d0
        enddo

        do i=2,N-1
            do j=2,M-1
                u(i,j) = 0.5d0*(psi(i,j+1)-psi(i,j-1))/dy
                v(i,j) = -0.5d0*(psi(i+1,j)-psi(i-1,j))/dx
            enddo
        enddo

        return
        end subroutine caluv


!!! compute temperature field
        subroutine calT(N,M,dt,dx,dy,u,v,T)
        implicit none
        integer :: i, j, N, M
        real(8) :: dx, dy, dt, dTx2, dTy2, dTx1, dTy1
        real(8) :: T(N,M), u(N,M), v(N,M), RT(N,M), Ti(N,M)

        call REST(N,M,dx,dy,dt,u,v,T,RT)
        do i=2,N-1
            do j=2,M-1
                Ti(i,j) = T(i,j)+dt*RT(i,j)
            enddo
        enddo
        call bcT(N,M,dx,dy,dt,Ti)

        call REST(N,M,dx,dy,dt,u,v,Ti,RT)
        do i=2,N-1
            do j=2,M-1
                Ti(i,j) = 0.75d0*T(i,j)+0.25d0*(Ti(i,j)+dt*RT(i,j))
            enddo
        enddo
        call bcT(N,M,dx,dy,dt,Ti)

        call REST(N,M,dx,dy,dt,u,v,Ti,RT)
        do i=2,N-1
            do j=2,M-1
                T(i,j) = 1.0d0/3.0d0*T(i,j)+2.0d0/3.0d0*(Ti(i,j)+dt*RT(i,j))
            enddo
        enddo
        call bcT(N,M,dx,dy,dt,T)

        return
        end subroutine calT

        subroutine REST(N,M,dx,dy,dt,u,v,T,RT)
        implicit none
        integer :: i, j, N, M
        real(8) :: dx, dy, dt, dTx2, dTy2, dTx1, dTy1
        real(8) :: T(N,M), u(N,M), v(N,M), RT(N,M)

       ! Interior points using FTCS Sheme
        do i=2,N-1
            do j=2,M-1
                dTx2 = (T(i+1,j)-2.0d0*T(i,j)+T(i-1,j))/dx/dx
                dTy2 = (T(i,j+1)-2.0d0*T(i,j)+T(i,j-1))/dy/dy
                dTx1 = (u(i+1,j)*T(i+1,j)-u(i-1,j)*T(i-1,j))/2.0/dx
                dTy1 = (v(i,j+1)*T(i,j+1)-v(i,j-1)*T(i,j-1))/2.0/dy
                RT(i,j) = dTx2+dTy2-dTx1-dTy1
            enddo
        enddo

        return
        end subroutine REST


        subroutine bcT(N,M,dx,dy,dt,T)
        implicit none
        integer :: i, j, N, M
        real(8) :: dx, dy, dt
        real(8) :: T(N,M)

       !Left and right side boundary(Dirichlet B.C.)
        do j=1,M
            T(1,j) = 1.0d0
            T(N,j) = 0.0d0
        enddo

       !Top and bottom side boundary(Neumann B.C.)
        do i=1,N
            T(i,1) = (4.0d0*T(i,2)-T(i,3))/3.0d0
            T(i,M) = (4.0d0*T(i,M-1)-T(i,M-2))/3.0d0
        enddo

        return
        end subroutine bcT



!!! check convergence
        subroutine convergence(N,M,dt,RVOR,Rpsi,error,itc)
        implicit none
        integer :: N, M, i, j, itc
        real(8) :: RVOR(N,M), Rpsi(N,M)
        real(8) :: dt, error, errvor, errpsi

        itc = itc+1
        errvor = 0.0d0
        errpsi = 0.0d0

        do i=1,N
            do j=1,M
                if(ABS(RVOR(i,j))*dt.GT.errvor) errvor = ABS(RVOR(i,j))*dt
                if(ABS(Rpsi(i,j)).GT.errpsi) errpsi = ABS(Rpsi(i,j))
            enddo
        enddo

        error = MAX(errvor,errpsi)
        if(itc.EQ.1) error = 100.0d0



        open(unit=01,file='error.dat',status='unknown',position='append')

        write(*,*) itc,' ',error
        if (MOD(itc,2000).EQ.0) then
            write(01,*) itc,' ',error
        endif

        close(01)


        return
        end subroutine convergence


!!! output data file
        subroutine output(N,M,X,Y,u,v,psi,T,k)
        implicit none
        integer :: N, M, i, j, k
        real(8) :: X(N), Y(M), u(N,M), v(N,M), psi(N,M),T(N,M)
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
                write(02,100) X(i), Y(j), u(i,j), v(i,j), psi(i,j), T(i,j)
            enddo
        enddo

100     format(2x,10(e12.6,'      '))
101     format('Title="Buoyancy Driven Cavity Flow"')
102     format('Variables=x,y,u,v,psi,T')
103     format('zone',1x,'i=',1x,i5,2x,'j=',1x,i5,1x,'f=point')

        close(02)

        return
        end subroutine output



