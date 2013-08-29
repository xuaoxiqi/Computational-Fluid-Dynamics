
!!!    This program sloves Lid Driven Cavity Flow problem using Vorticity-Streamfunction Method
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
        integer, parameter :: nx=129, ny=129
        integer :: i, j, itc, itc_max
        real(8) :: dx, dy, Re, dt, eps, error
        real(8) :: X(nx), Y(ny), u(nx,ny), v(nx,ny), vor(nx,ny), RVOR(nx,ny), psi(nx,ny), Rpsi(nx,ny)

!!! input initial data
        Re = 1000.0d0
        dx = 1.0d0/(nx-1)
        dy = 1.0d0/(ny-1)
        dt = 1e-4
        eps = 1e-5
        itc = 0
        itc_max = 5*1e5
        error = 100.0d0

!!! set up initial flow field
        call initial(nx,ny,dx,dy,X,Y,u,v,psi,vor,RVOR,Rpsi)

        do while((error.GT.eps).AND.(itc.LT.itc_max))

!!! solve Vorticity Equation
            call solvor(nx,ny,dx,dy,Re,dt,u,v,vor,RVOR)

!!! solve Streamfunction Equation
            call solpsi(nx,ny,dx,dy,vor,psi,Rpsi)

!!! updates the values of sream function at boundary points
            call bcpsi(nx,ny,dy,psi)

!!! updates the boundary condition for vorticity
            call bcvor(nx,ny,dx,dy,vor,psi)

!!! compute the velocity components u and v
            call caluv(nx,ny,dx,dy,psi,u,v)

!!! check convergence
            call convergence(nx,ny,dt,RVOR,Rpsi,error,itc)

!!! output preliminary results
            if (MOD(itc,20000).EQ.0) then
                call output(nx,ny,X,Y,u,v,psi,itc)
            endif

        enddo

!!! output data file
        call output(nx,ny,X,Y,u,v,psi,itc)

        write(*,*)
        write(*,*) '************************************************************'
        write(*,*) 'This program sloves Lid Driven Cavity Flow problem'
        write(*,*) 'using Vorticity-Streamfunction Methods'
        write(*,*) 'nx =',nx,',       ny =',ny
        write(*,*) 'Re=',Re
        write(*,*) 'dt =', dt
        write(*,*) 'eps =',eps
        write(*,*) 'itc =',itc
        write(*,*) 'Developing time=',dt*itc,'s'
        write(*,*) '************************************************************'
        write(*,*)

        stop
        end program main


!!! set up initial flow field
        subroutine initial(nx,ny,dx,dy,X,Y,u,v,psi,vor,RVOR,Rpsi)
        implicit none
        integer :: i, j, nx, ny
        real(8) :: dx, dy
        real(8) :: X(nx), Y(ny), u(nx,ny), v(nx,ny), psi(nx,ny), vor(nx,ny), RVOR(nx,ny), Rpsi(nx,ny)

        do i=1,nx
            X(i) = (i-1)*dx
        enddo
        do j=1,ny
            Y(j) = (j-1)*dy
        enddo

        u = 0.0d0
        v = 0.0d0
        psi = 0.0d0
        vor = 0.0d0
        RVOR = 0.0d0
        Rpsi = 0.0d0

        do i=1,nx
            u(i,ny) = 1.0d0
        enddo

        return
        end subroutine initial


!!! solve Vorticity Equation
        subroutine solvor(nx,ny,dx,dy,Re,dt,u,v,vor,RVOR)
        implicit none
        integer :: i, j, nx, ny
        real(8) :: dx, dy, dt, Re, dvorx2, dvory2, dvorx1, dvory1
        real(8) :: vor(nx,ny), u(nx,ny), v(nx,ny), RVOR(nx,ny)

       ! FTCS Sheme
        do i=2,nx-1
            do j=2,ny-1
                dvorx2 = (vor(i+1,j)-2.0d0*vor(i,j)+vor(i-1,j))/dx/dx
                dvory2 = (vor(i,j+1)-2.0d0*vor(i,j)+vor(i,j-1))/dy/dy
                dvorx1 = (u(i+1,j)*vor(i+1,j)-u(i-1,j)*vor(i-1,j))/2.0d0/dx
                dvory1 = (v(i,j+1)*vor(i,j+1)-v(i,j-1)*vor(i,j-1))/2.0d0/dy
                RVOR(i,j) = (dvorx2+dvory2)/Re-dvorx1-dvory1
                vor(i,j) = vor(i,j)+dt*RVOR(i,j)
            enddo
        enddo

        return
        end subroutine solvor


!!! solve Streamfunction Equation
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
                S(i,j)=vor(i,j)-(psi(i+1,j)-2.0d0*psi(i,j)+psi(i-1,j))/dx/dx-(psi(i,j+1)-2.0d0*psi(i,j)+psi(i,j-1))/dy/dy
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

        alpha = 1.5d0       !alpha is ralaxtion factor
        do i=3,nx-2
            do j=3,ny-2
                Rpsi(i,j)=(S(i,j)-aw*Rpsi(i-1,j)-as*Rpsi(i,j-1))/ap
                psi(i,j) = psi(i,j)+alpha*Rpsi(i,j)
            enddo
        enddo

        return
        end subroutine solpsi


!!! updates the values of Sreamfunction at boundary points
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
            psi(i,ny-1) = 0.25d0*(psi(i,ny-2)-2.0d0*dy)
        enddo

        return
        end subroutine bcpsi


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
        do i=2,nx-1
            vor(i,1) = 3.0d0*psi(i,2)/dy/dy-0.5d0*vor(i,2)
            vor(i,ny) = 3.0d0*(psi(i,ny-1)+dy)/dy/dy-0.5d0*vor(i,ny-1)
        enddo

        return
        end subroutine bcvor


!!! compute the velocity components u and v
        subroutine caluv(nx,ny,dx,dy,psi,u,v)
        implicit none
        integer :: i, j, nx, ny
        real(8) :: dx, dy
        real(8) :: psi(nx,ny), u(nx,ny), v(nx,ny)

       !physical boundary condition
        do j=1,ny
            u(1,j) = 0.0d0
            u(nx,j) = 0.0d0
            v(1,j) = 0.0d0
            v(nx,j) = 0.0d0
        enddo
        do i=2,nx-1
            u(i,1) = 0.0d0
            v(i,1) = 0.0d0
            u(i,ny) = 1.0d0
            v(i,ny) = 0.0d0
        enddo

        do i=2,nx-1
            do j=2,ny-1
                u(i,j) = 0.5d0*(psi(i,j+1)-psi(i,j-1))/dy
                v(i,j) = -0.5d0*(psi(i+1,j)-psi(i-1,j))/dx
            enddo
        enddo

        return
        end subroutine caluv


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

!        write(*,*) itc,' ',error

!!!        open(unit=01,file='error.dat',status='unknown',position='append')
!!!        if (MOD(itc,2000).EQ.0) then
!!!            write(01,*) itc,' ',error
!!!        endif
!!!        close(01)

        return
        end subroutine convergence



!!! output data file
    subroutine output(nx,ny,X,Y,uc,vc,psi,itc)
    implicit none
    integer :: nx, ny, i, j, itc
    real(8) :: X(nx), Y(ny), uc(nx,ny), vc(nx,ny), psi(nx,ny)
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
