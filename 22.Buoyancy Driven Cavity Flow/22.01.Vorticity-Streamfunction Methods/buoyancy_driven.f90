

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
!!!    Phi(i,j)-------------------stream function
!!!    vor(i,j)-------------------vorticity function
!!!    RPhi(i,j)------------------Phi^{n+1}_{i,j} - Phi^{n}_{i,j}
!!!    RVOR(i,j)------------------(vor^{n+1}_{i,j}-vor^{n}_{i,j})/dt
!!!    T(i,j))--------------------Temperature field
!!!    RT(i,j)--------------------(T^{n+1}_{i,j}-T^{n}_{i,j})/dt

       program main
       implicit none
       integer, parameter :: N=201, M=201
       integer :: i, j, itc
       real :: dx, dy, Pr, Ra, dt, eps, errvor, errphi
       real :: X(N), Y(M), u(N,M), v(N,M), vor(N,M), RVOR(N,M), Phi(N,M), RPhi(N,M), T(N,M), RT(N,M)

!!! input initial data
       Pr = 0.71
       Ra = 1e3
       dx = 1./(N-1)
       dy = 1./(M-1)
       dt = 1e-6
       eps = 1e-4
       itc = 0

       write(*,*) 'This program sloves Buoyancy Driven Cavity Flow problem using Vorticity-Streamfunction Methods'
!!! set up initial flow field
       call initial(N,M,dx,dy,X,Y,u,v,Phi,vor,RVOR,RPhi,T,RT)

       do
!!! solve vorticity equation
              call solvor(N,M,dx,dy,Pr,Ra,dt,u,v,vor,RVOR,T)

!!! solve stream function equation
              call solphi(N,M,dx,dy,vor,Phi,RPhi)

!!! updates the values of sream function at boundary points
              call bcphi(N,M,dy,Phi)

!!! updates the boundary condition for vorticity
              call bcvor(N,M,dx,dy,vor,Phi)

!!! compute the velocity components u and v
              call caluv(N,M,dx,dy,Phi,u,v)

!!! compute the temperature field
              call calt(N,M,dt,dx,dy,u,v,T,RT)

              errvor = 0.0
              errphi = 0.0
              do i=1,N
                     do j=1,M
                            if(ABS(RVOR(i,j))*dt.GT.errvor) errvor = ABS(RVOR(i,j))*dt
                            if(ABS(RPhi(i,j)).GT.errphi) errphi = ABS(RPhi(i,j))
                     enddo
              enddo
              !print*, 'errvor=',errvor
              !print*, 'errphi=',errphi

              if((MAX(errvor,errphi).LT.eps).AND.(itc.GT.1).OR.(itc.EQ.5*1e5)) then
                     write(*,*) 'Results meet convergence criteria!'
                     exit
              else
                     itc = itc+1
                     !print*,'Iterative times=',itc
                     cycle
              endif

       enddo

       write(*,*) 'eps=',eps
       write(*,*) 'Pr=',Pr
       write(*,*) 'Ra=',Ra
       write(*,*) 'dt =', dt
       write(*,*) 'Iterative times=',itc
       write(*,*) 'Developing time=',dt*itc

       open(unit=02,file='./cavity.dat',status='unknown')
       write(02,101)
       write(02,102)
       write(02,103) N, M
       do j=1,M
              do i = 1,N
                     write(02,100) X(i), Y(j), u(i,j), v(i,j), Phi(i,j), VOR(i,j), T(i,j)
              enddo
       enddo

100    format(2x,10(e12.6,'      '))
101    format('Title="Driven Cavity Flow"')
102    format('Variables=x,y,u,v,Phi,VOR,T')
103    format('zone',1x'i=',1x,i5,2x,'j=',1x,i5,1x,'f=point')

       close(02)
       write(*,*) 'Data export to ./cavity.dat file!'

       stop
       end program main



       subroutine initial(N,M,dx,dy,X,Y,u,v,Phi,vor,RVOR,RPhi,T,RT)
       implicit none
       integer :: i, j, N, M
       real :: dx, dy
       real :: X(N), Y(M), u(N,M), v(N,M), Phi(N,M), vor(N,M), RVOR(N,M), RPhi(N,M), T(N,M), RT(N,M)

       do i=1,N
              X(i) = (i-1)*dx
       enddo
       do j=1,M
              Y(j) = (j-1)*dy
       enddo

       do i=1,N
              do j=1,M
                     u(i,j) = 0.0
                     v(i,j) = 0.0!!!    u(i,j), v(i,j)------------velocity function
                     Phi(i,j) = 0.0!!!    Phi(i,j)--------------------stream function
                     RPhi(i,j) = 0.0!!!    RPhi(i,j)--------------Phi^{n+1}_{i,j} - Phi^{n}_{i,j}
                     vor(i,j) = 0.0!!!    vor(i,j)------------------vorticity function
                     RVOR(i,j) = 0.0!!!    RVOR(i,j)-------------(vor^{n+1}_{i,j}-vor^{n}_{i,j})/dt
                     T(i,j) = 0.0!!!      T(i,j))--------------Temperature field
                     RT(i,j) = 0.0!!!


              enddo
       enddo

       do j =1,M
              T(1,j) = 1.0  !Left side hot wall
       enddo

       return
       end subroutine initial



       subroutine solvor(N,M,dx,dy,Pr,Ra,dt,u,v,vor,RVOR,T)
       implicit none
       integer :: i, j, N, M
       real :: dx, dy, dt, Pr, Ra, dvorx2, dvory2, dvorx1, dvory1
       real :: vor(N,M), u(N,M), v(N,M), RVOR(N,M), T(N,M)

       ! FTCS Sheme
       do i=2,N-1
              do j=2,M-1
                     dvorx2 = (vor(i+1,j)-2*vor(i,j)+vor(i-1,j))/dx/dx
                     dvory2 = (vor(i,j+1)-2*vor(i,j)+vor(i,j-1))/dy/dy
                     dvorx1 = (u(i+1,j)*vor(i+1,j)-u(i-1,j)*vor(i-1,j))/2/dx
                     dvory1 = (v(i,j+1)*vor(i,j+1)-v(i,j-1)*vor(i,j-1))/2/dy
                     RVOR(i,j) = (dvorx2+dvory2)*Pr-dvorx1-dvory1+Ra*Pr*(T(i+1,j)-T(i-1,j))/2/dx
                     vor(i,j) = vor(i,j)+dt*RVOR(i,j)
              enddo
       enddo

       return
       end subroutine solvor



       subroutine solphi(N,M,dx,dy,vor,Phi,RPhi)
       implicit none
       integer :: i, j ,N, M
       real :: alpha, dx, dy, aw, as, ap
       real :: vor(N,M), Phi(N,M), RPhi(N,M), S(N,M)

       aw = 1.0/dx/dx
       as = 1.0/dy/dy
       ap = -2*(as+aw)

       do i=3,N-2
              do j=3,M-2
                     S(i,j)=vor(i,j)-(Phi(i+1,j)-2*Phi(i,j)+Phi(i-1,j))/dx/dx-(Phi(i,j+1)-2*Phi(i,j)+Phi(i,j-1))/dy/dy
              enddo
       enddo

       do j=1,M
              RPhi(1,j) = 0.0
              RPhi(2,j) = 0.0
              RPhi(N,j) = 0.0
              RPhi(N-1,j) = 0.0
       enddo
       do i=1,N
              RPhi(i,1) = 0.0
              RPhi(i,2) = 0.0
              RPhi(i,M) = 0.0
              RPhi(i,M-1) = 0.0
       enddo

       alpha = 0.7      !alpha is ralaxtion factor
       do i=3,N-2
              do j=3,M-2
                     RPhi(i,j)=(S(i,j)-aw*RPhi(i-1,j)-as*RPhi(i,j-1))/ap
                     Phi(i,j) = Phi(i,j)+alpha*RPhi(i,j)
              enddo
       enddo

       return
       end subroutine solphi



       subroutine bcphi(N,M,dy,Phi)
       implicit none
       integer :: i, j, N, M
       real :: dy
       real :: Phi(N,M)

       do j=2,M-1
              Phi(2,j) = 0.25*Phi(3,j)
              Phi(N-1,j) = 0.25*Phi(N-2,j)
       enddo
       do i=2,N-1
              Phi(i,2) = 0.25*Phi(i,3)
              Phi(i,M-1) = 0.25*(Phi(i,M-2)-2.0*dy)
       enddo

       return
       end subroutine bcphi



       subroutine bcvor(N,M,dx,dy,vor,Phi)
       implicit none
       integer :: i, j, N, M
       real :: dx, dy
       real :: vor(N,M), Phi(N,M)

       ! 2nd order approximation
       do j=1,M
              vor(1,j) = 3.0*Phi(2,j)/dx/dx-0.5*vor(2,j)
              vor(N,j) = 3.0*Phi(N-1,j)/dx/dx-0.5*vor(N-1,j)
       enddo
       do i=2,N-1
              vor(i,1) = 3.0*Phi(i,2)/dy/dy-0.5*vor(i,2)
              vor(i,M) = 3.0*(Phi(i,M-1)+dy)/dy/dy-0.5*vor(i,M-1)
       enddo

       return
       end subroutine bcvor



       subroutine caluv(N,M,dx,dy,Phi,u,v)
       implicit none
       integer :: i, j, N, M
       real :: dx, dy
       real :: Phi(N,M), u(N,M), v(N,M)

       !physical boundary condition
       do j=1,M
              u(1,j) = 0
              u(N,j) = 0
              v(1,j) = 0
              v(N,j) = 0
       enddo
       do i=2,N-1
              u(i,1) = 0
              v(i,1) = 0
              u(i,M) = 1
              v(i,M) = 0
       enddo

       do i=2,N-1
              do j=2,M-1
                     u(i,j) = 0.5*(Phi(i,j+1)-Phi(i,j-1))/dy
                     v(i,j) = -0.5*(Phi(i+1,j)-Phi(i-1,j))/dx
              enddo
       enddo

       return
       end subroutine caluv



       subroutine calt(N,M,dt,dx,dy,u,v,T,RT)
       implicit none
       integer :: i, j, N, M
       real :: dx, dy, dt, dTx2, dTy2, dTx1, dTy1
       real :: T(N,M), u(N,M), v(N,M), RT(N,M)

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
              T(1,j) = 1.0
              T(N,j) = 0.0
       enddo

       !Top and bottom side boundary(Neumann B.C.)
       do i=1,N
              T(i,1) = (4*T(i,2)-T(i,3))/3.0
              T(i,M) = (4*T(i,M-1)-T(i,M-2))/3.0
       enddo


       end subroutine calt

