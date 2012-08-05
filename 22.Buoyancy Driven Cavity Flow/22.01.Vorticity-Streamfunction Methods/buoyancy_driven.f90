

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

!!!    u(i,j), v(i,j)-------------velocity function
!!!    psi(i,j)-------------------stream function
!!!    vor(i,j)-------------------vorticity function
!!!    Rpsi(i,j)------------------psi^{n+1}_{i,j} - psi^{n}_{i,j}
!!!    RVOR(i,j)------------------(vor^{n+1}_{i,j}-vor^{n}_{i,j})/dt
!!!    T(i,j))--------------------Temperature field
!!!    RT(i,j)--------------------(T^{n+1}_{i,j}-T^{n}_{i,j})/dt

        program main
        implicit none
        integer, parameter :: N=101, M=101
        integer :: i, j, itc, itc_max
        real(8) :: dx, dy, Pr, Ra, dt, eps, error
        real(8) :: X(N), Y(M), u(N,M), v(N,M), vor(N,M), RVOR(N,M), psi(N,M), Rpsi(N,M), T(N,M)
        real(8) :: psi_mid, u_max, v_max, x_loc, y_loc

!!! input initial data
        Pr = 0.71d0
        Ra = 1e4
        dx = 1.0d0/(N-1)
        dy = 1.0d0/(M-1)
        dt = 1e-7
        eps = 1e-3
        itc = 0
        itc_max = 5*1e6
        error = 100.0d0

        write(*,*) 'Start:'
!!! set up initial flow field
        call initial(N,M,dx,dy,X,Y,u,v,psi,vor,RVOR,Rpsi,T)

        do while((error.GT.eps).AND.(itc.LT.itc_max))
!!! solve vorticity equation
            call solvor(N,M,dx,dy,Pr,Ra,dt,u,v,vor,RVOR,T)

!!! solve Streamfunction equation
            call solpsi(N,M,dx,dy,vor,psi,Rpsi)

!!! updates the values of sream function at boundary points
            call bcpsi(N,M,dy,psi)

!!! updates the boundary condition for vorticity
            call bcvor(N,M,dx,dy,vor,psi)

!!! compute velocity components u and v
            call caluv(N,M,dx,dy,psi,u,v)

!!! compute temperature field
            call calt(N,M,dt,dx,dy,u,v,T)

!!! check convergence
            call convergence(N,M,dt,RVOR,Rpsi,error,itc)

        enddo

!!! validate results with reference
        call check(N,M,dx,dy,psi,u,v,psi_mid,u_max,v_max,x_loc,y_loc)

!!! output data file
        call output(N,M,X,Y,u,v,psi,VOR,T)

        write(*,*)
        write(*,*) '************************************************************'
        write(*,*) 'This program sloves Buoyancy Driven Cavity Flow problem'
        write(*,*) 'using Vorticity-Streamfunction Methods'
        write(*,*) 'N =',N,',       M =',M
        write(*,*) 'Pr=',Pr
        write(*,*) 'Ra=',Ra
        write(*,*) 'dt =', dt
        write(*,*) 'eps =',eps
        write(*,*) 'itc =',itc
        write(*,*) 'Developing time=',dt*itc,'s'
        write(*,*)
        write(*,*) 'psi_mid =',psi_mid
        write(*,*) 'u_max =',u_max,'at y =',y_loc
        write(*,*) 'v_max =',v_max,'at x =',x_loc
        write(*,*) '************************************************************'
        write(*,*)

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

        do i=1,N
            do j=1,M
                u(i,j) = 0.0d0
                v(i,j) = 0.0d0!!!    u(i,j), v(i,j)------------velocity function
                psi(i,j) = 0.0d0!!!    psi(i,j)--------------------stream function
                Rpsi(i,j) = 0.0d0!!!    Rpsi(i,j)--------------psi^{n+1}_{i,j} - psi^{n}_{i,j}
                vor(i,j) = 0.0d0!!!    vor(i,j)------------------vorticity function
                RVOR(i,j) = 0.0d0!!!    RVOR(i,j)-------------(vor^{n+1}_{i,j}-vor^{n}_{i,j})/dt
                T(i,j) = 0.0d0!!!      T(i,j))--------------Temperature field
            enddo
        enddo

        do j =1,M
            T(1,j) = 1.0d0  !Left side hot wall
        enddo

        return
        end subroutine initial


!!! solve vorticity equation
        subroutine solvor(N,M,dx,dy,Pr,Ra,dt,u,v,vor,RVOR,T)
        implicit none
        integer :: i, j, N, M
        real(8) :: dx, dy, dt, Pr, Ra, dvorx2, dvory2, dvorx1, dvory1
        real(8) :: vor(N,M), u(N,M), v(N,M), RVOR(N,M), T(N,M)

        ! FTCS Sheme
        do i=2,N-1
            do j=2,M-1
                dvorx2 = (vor(i+1,j)-2.0d0*vor(i,j)+vor(i-1,j))/dx/dx
                dvory2 = (vor(i,j+1)-2.0d0*vor(i,j)+vor(i,j-1))/dy/dy
                dvorx1 = (u(i+1,j)*vor(i+1,j)-u(i-1,j)*vor(i-1,j))/2.0d0/dx
                dvory1 = (v(i,j+1)*vor(i,j+1)-v(i,j-1)*vor(i,j-1))/2.0d0/dy
                RVOR(i,j) = (dvorx2+dvory2)*Pr-dvorx1-dvory1+Ra*Pr*(T(i+1,j)-T(i-1,j))/2.0d0/dx
                vor(i,j) = vor(i,j)+dt*RVOR(i,j)
            enddo
        enddo

        return
        end subroutine solvor


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
                u(i,j) = -0.5d0*(psi(i,j+1)-psi(i,j-1))/dy
                v(i,j) = 0.5d0*(psi(i+1,j)-psi(i-1,j))/dx
            enddo
        enddo

        return
        end subroutine caluv


!!! compute temperature field
        subroutine calt(N,M,dt,dx,dy,u,v,T)
        implicit none
        integer :: i, j, N, M
        real(8) :: dx, dy, dt, dTx2, dTy2, dTx1, dTy1
        real(8) :: T(N,M), u(N,M), v(N,M), RT(N,M)

       ! Interior points using FTCS Sheme
        do i=2,N-1
            do j=2,M-1
                dTx2 = (T(i+1,j)-2*T(i,j)+T(i-1,j))/dx/dx
                dTy2 = (T(i,j+1)-2*T(i,j)+T(i,j-1))/dy/dy
                dTx1 = (u(i+1,j)*T(i+1,j)-u(i-1,j)*T(i-1,j))/2/dx
                dTy1 = (v(i,j+1)*T(i,j+1)-v(i,j-1)*T(i,j-1))/2/dy
                RT(i,j) = dTx2+dTy2+dTx1+dTy1
                T(i,j) = T(i,j)+dt*RT(i,j)
            enddo
        enddo

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
        end subroutine calt

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

        write(*,*) 'itc=',itc,'    |    error=',error

        return
        end subroutine convergence


!!! validate results with reference
        subroutine check(N,M,dx,dy,psi,u,v,psi_mid,u_max,v_max,x_loc,y_loc)
        implicit none
        integer :: N, M, mid_x, mid_y, i, j, temp_x, temp_y
        real(8) :: dx, dy, psi_mid, u_max, v_max, x_loc, y_loc
        real(8) :: psi(N,M), u(N,M), v(N,M)

        mid_x = INT(N/2)
        mid_y = INT(M/2)
        psi_mid = psi(mid_x,mid_y)

        u_max = 0.0d0
        v_max = 0.0d0
        temp_x = 0
        temp_y = 0
        do j=1,M
            if(u(mid_x,j).GT.u_max) then
                u_max = u(mid_x,j)
                temp_y = j
            endif
        enddo
        y_loc = (temp_y-1)*dy

        do i=1,N
            if(v(i,mid_y).GT.v_max) then
                v_max = v(i,mid_y)
                temp_x = i
            endif
        enddo
        x_loc = (temp_x-1)*dx

        return
        end subroutine check

!!! output data file
        subroutine output(N,M,X,Y,u,v,psi,VOR,T)
        implicit none
        integer :: N, M, i, j
        real(8) :: X(N), Y(M), u(N,M), v(N,M), psi(N,M), VOR(N,M), T(N,M)

        open(unit=02,file='./cavity.dat',status='unknown')
        write(02,101)
        write(02,102)
        write(02,103) N, M
        do j=1,M
            do i = 1,N
                write(02,100) X(i), Y(j), u(i,j), v(i,j), psi(i,j), VOR(i,j), T(i,j)
            enddo
        enddo

100     format(2x,10(e12.6,'      '))
101     format('Title="Lid Driven Cavity Flow(Vorticity-Streamfunction Methods)"')
102     format('Variables=x,y,u,v,psi,VOR,T')
103     format('zone',1x,'i=',1x,i5,2x,'j=',1x,i5,1x,'f=point')

        close(02)
        write(*,*) 'Data export to ./cavity.dat file!'

        return
        end subroutine output



