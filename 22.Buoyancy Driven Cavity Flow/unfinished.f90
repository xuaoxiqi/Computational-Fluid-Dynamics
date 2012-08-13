
!!!    This program sloves Buoyancy Driven Cavity Flow problem using Artificial Compressibility Methods
!!!    Solve Momentum Equation with QUICK Scheme
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
        integer, parameter :: N=101,M=101
        integer :: itc, itc_max, i, j, k
        real(8) :: u(N,M+1),v(N+1,M),p(N+1,M+1),T(N+1,M+1), psi(N,M), X(N), Y(M)
        real(8) :: un(N,M+1),vn(N+1,M),pn(N+1,M+1),Tn(N+1,M+1),uc(N,M),vc(N,M),pc(N,M),Tc(N,M)
        real(8) :: c, c2, Ra, Pr, dt, dx, dy, eps, error
        real(8) :: psi_mid, u_max, v_max, x_loc, y_loc

!!! input initial data
        c = 1.5d0
        c2 = c*c   ! c2 = 2.25d0
        Ra = 1e4
        Pr = 0.71d0
        dx = 1.0d0/float(N-1)
        dy = 1.0d0/float(M-1)
        dt = 2*1e-5
        eps = 1e-6
        itc_max = 2*1e6
        itc = 0
        error = 100.00d0
        k = 0

        write(*,*) 'Start:'
!!! set up initial flow field
        call initial(N,M,dx,dy,X,Y,u,v,p,T,psi)

        do while((error.GT.eps).AND.(itc.LT.itc_max))

            error=0.0d0

!!! Solve Momentum Equation with QUICK Scheme
            call quick(N,M,dx,dy,dt,Ra,Pr,u,v,p,T,un,vn)

!!! Solve Continuity Equation
            call calpn(N,M,dx,dy,dt,c2,p,un,vn,pn)

!!! Solve Energy Equation
            call caltn(N,M,dx,dy,dt,T,un,vn,Tn)

!!! check convergence
            call convergence(N,M,dt,c2,error,u,v,p,T,un,vn,pn,Tn,itc)

!!! output preliminary results
            if (MOD(itc,10000).EQ.0) then
                call caluvpt(N,M,u,v,p,T,uc,vc,pc,Tc)
                call calpsi(N,M,dx,dy,uc,vc,psi)
                k = k+1
                call output(N,M,X,Y,uc,vc,psi,pc,Tc,k)
            endif
            
        enddo

!!! compute velocity components u, v and pressure p
        call caluvpt(N,M,u,v,p,T,uc,vc,pc,Tc)

!!! compute Streamfunction
        call calpsi(N,M,dx,dy,uc,vc,psi)

!!! validate results with reference
        call check(N,M,dx,dy,psi,u,v,psi_mid,u_max,v_max,x_loc,y_loc)

!!! output data file
        call output(N,M,X,Y,uc,vc,psi,pc,Tc,k)

        write(*,*)
        write(*,*) '************************************************************'
        write(*,*) 'This program sloves Buoyancy Driven Cavity Flow problem'
        write(*,*) 'using Artificial Compressibility Methods'
        write(*,*) 'N =',N,',       M =',M
        write(*,*) 'Ra =',Ra
        write(*,*) 'Rr =',Pr
        write(*,*) 'dt =',dt
        write(*,*) 'c (Artificial Compressibility coefficient) =',c
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
        subroutine initial(N,M,dx,dy,X,Y,u,v,p,T,psi)
        implicit none
        integer :: N, M, i, j
        real(8) :: dx, dy
        real(8) :: u(N,M+1), v(N+1,M), p(N+1,M+1), T(N+1,M+1), psi(N,M), X(N), Y(M)

        do i=1,N
            X(i) = (i-1)*dx
        enddo
        do j=1,M
            Y(j) = (j-1)*dy
        enddo
        do i=1,N+1
            do j=1,M+1
                p(i,j) = 1.0d0
                T(i,j) = 0.0d0
                if(i.EQ.1) T(i,j) = 4.0d0/3.0d0
                if(i.EQ.2) T(i,j) = 2.0d0/3.0d0
            enddo
        enddo
        do i=1,N
            do j=1,M+1
                u(i,j) = 0.0
            enddo
        enddo
        do i=1,N+1
            do j=1,M
                v(i,j) = 0.0d0
            enddo
        enddo

        do i=1,N
            do j=1,M
                psi(i,j) = 0.0d0
            enddo
        enddo

        return
        end subroutine initial


!!! Solve Momentum Equation with QUICK Scheme
        subroutine quick(N,M,dx,dy,dt,Ra,Pr,u,v,p,T,un,vn)
        implicit none
        integer :: N, M, i, j
        real(8) :: u(N,M+1),v(N+1,M),p(N+1,M+1),T(N+1,M+1),un(N,M+1),vn(N+1,M)
        real(8) :: Ra, Pr, dx, dy, dt
        real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
        real(8) :: alpha

!!!!!!!!!!!!!!!!!!!!!compute x-direction velocity component un!!!!!!!!!!!!!!!!!!!!!
        do i=3,N-2
            do j=3,M-1

                fw = 0.5d0*(u(i-1,j)+u(i,j))*dy
                fe = 0.5d0*(u(i,j)+u(i+1,j))*dy
                fs = 0.5d0*(v(i,j-1)+v(i+1,j-1))*dx
                fn = 0.5d0*(v(i,j)+v(i+1,j))*dx
                df = fe-fw+fn-fs

                !!! common coefficient in 3rd-order upwind QUICK Scheme
                aw = Pr+0.75d0*alpha(fw)*fw+0.125d0*alpha(fe)*fe+0.375d0*(1.0d0-alpha(fw))*fw
                aww = -0.125d0*alpha(fw)*fw
                ae = Pr-0.375d0*alpha(fe)*fe-0.75d0*(1.0-alpha(fe))*fe-0.125d0*(1.0d0-alpha(fw))*fw
                aee = 0.125d0*(1.0d0-alpha(fe))*fe
                as = Pr+0.75d0*alpha(fs)*fs+0.125d0*alpha(fn)*fn+0.375d0*(1.0d0-alpha(fs))*fs
                ass = -0.125d0*alpha(fs)*fs
                an = Pr-0.375d0*alpha(fn)*fn-0.75d0*(1.0-alpha(fn))*fn-0.125d0*(1.0d0-alpha(fs))*fs
                ann = 0.125d0*(1.0d0-alpha(fn))*fn
                ap = aw+ae+as+an+aww+aee+ass+ann+df

                un(i,j) = u(i,j) + dt/dx/dy*( -ap*u(i,j)+aw*u(i-1,j)+ae*u(i+1,j)+aww*u(i-2,j)+aee*u(i+2,j)&
                                +as*u(i,j-1)+an*u(i,j+1)+ass*u(i,j-2)+ann*u(i,j+2) ) - dt*(p(i+1,j)-p(i,j))/dx

            enddo
        enddo

        !!! compute interior region boundary with 1st-order upwind discrete scheme
        j=2
        do i=3,N-2
            call upbound_u(N,M,dx,dy,dt,Ra,Pr,u,v,p,un,i,j)
        enddo
        j=M
        do i=3,N-2
            call upbound_u(N,M,dx,dy,dt,Ra,Pr,u,v,p,un,i,j)
        enddo
        i=2
        do j=2,M
            call upbound_u(N,M,dx,dy,dt,Ra,Pr,u,v,p,un,i,j)
        enddo
        i=N-1
        do j=2,M
            call upbound_u(N,M,dx,dy,dt,Ra,Pr,u,v,p,un,i,j)
        enddo

        !!! compute exterior region boundary with physical boundary condition
        do i=2,N-1
            un(i,1) = -un(i,2)
            un(i,M+1) = -un(i,M)
        enddo
        do j=1,M+1
            un(1,j) = 0.0d0
            un(N,j) = 0.0d0
        enddo
!!!!!!!!!!!!!!!!!!!!!compute x-direction velocity component un!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!compute y-direction velocity component vn!!!!!!!!!!!!!!!!!!!!!
        do i=3,N-1
            do j=3,M-2

                fw = 0.5d0*(u(i-1,j)+u(i-1,j+1))*dy
                fe = 0.5d0*(u(i,j)+u(i,j+1))*dy
                fs = 0.5d0*(v(i,j-1)+v(i,j))*dx
                fn = 0.5d0*(v(i,j)+v(i,j+1))*dx
                df = fe-fw+fn-fs

                !!! common coefficient in 3rd-order upwind QUICK Scheme
                aw = Pr+0.75d0*alpha(fw)*fw+0.125d0*alpha(fe)*fe+0.375d0*(1.0d0-alpha(fw))*fw
                aww = -0.125d0*alpha(fw)*fw
                ae = Pr-0.375d0*alpha(fe)*fe-0.75d0*(1.0-alpha(fe))*fe-0.125d0*(1.0d0-alpha(fw))*fw
                aee = 0.125d0*(1.0d0-alpha(fe))*fe
                as = Pr+0.75d0*alpha(fs)*fs+0.125d0*alpha(fn)*fn+0.375d0*(1.0d0-alpha(fs))*fs
                ass = -0.125d0*alpha(fs)*fs
                an = Pr-0.375d0*alpha(fn)*fn-0.75d0*(1.0-alpha(fn))*fn-0.125d0*(1.0d0-alpha(fs))*fs
                ann = 0.125d0*(1.0d0-alpha(fn))*fn
                ap = aw+ae+as+an+aww+aee+ass+ann+df

                vn(i,j) = v(i,j) + dt/dx/dy*( -ap*v(i,j)+aw*v(i-1,j)+ae*v(i+1,j)+aww*v(i-2,j)+aee*v(i+2,j)   &
                                +as*v(i,j-1)+an*v(i,j+1)+ass*v(i,j-2)+ann*v(i,j+2) ) - dt*(p(i,j+1)-p(i,j))/dy &
                                + dt*Ra*Pr*(T(i,j)+T(i,j+1))/2.0d0

            enddo
        enddo

        !!! compute interior region boundary with 1st-order upwind discrete scheme
        j=2
        do i=3,N-1
            call upbound_v(N,M,dx,dy,dt,Ra,Pr,u,v,p,T,vn,i,j)
        enddo
        j=M-1
        do i=3,N-1
            call upbound_v(N,M,dx,dy,dt,Ra,Pr,u,v,p,T,vn,i,j)
        enddo
        i=2
        do j=2,M-1
            call upbound_v(N,M,dx,dy,dt,Ra,Pr,u,v,p,T,vn,i,j)
        enddo
        i=N
        do j=2,M-1
            call upbound_v(N,M,dx,dy,dt,Ra,Pr,u,v,p,T,vn,i,j)
        enddo

       !!! compute exterior region boundary with physical boundary condition
        do i=2,N
            vn(i,1) = 0.0d0
            vn(i,M) = 0.0d0
        enddo
        do j=1,M
            vn(1,j) = -vn(2,j)
            vn(N+1,j) = -vn(N,j)
        enddo
!!!!!!!!!!!!!!!!!!!!!compute y-direction velocity component vn!!!!!!!!!!!!!!!!!!!!!

        return
        end subroutine quick


!!! if(f_k.GT.0) then alpha_k = 1   (k=w,e,s,n)
!!! if(f_k.LT.0) then alpha_k = 0   (k=w,e,s,n)
        function alpha(x)
        implicit none
        real(8) :: alpha, x

        if(x.GT.0.0d0) then
            alpha = 1.0d0
        elseif(x.LT.0.0d0) then
            alpha = 0.0d0
        endif

        return
        end function alpha


!!! compute interior region boundary with 1st-order upwind discrete scheme-->un
        subroutine upbound_u(N,M,dx,dy,dt,Ra,Pr,u,v,p,un,i,j)
        implicit none
        integer :: N, M, i, j
        real(8) :: dx, dy, dt, Ra, Pr
        real(8) :: u(N,M+1),v(N+1,M),p(N+1,M+1),un(N,M+1)
        real(8) :: aw, ae, as, an, df, ap

        aw = Pr+MAX(0.5d0*(u(i-1,j)+u(i,j))*dy,0.0d0)
        ae = Pr+MAX(0.0d0,-0.5d0*(u(i,j)+u(i+1,j))*dy)
        as = Pr+MAX(0.5d0*(v(i,j-1)+v(i+1,j-1))*dx,0.0d0)
        an = Pr+MAX(0.0d0,-0.5d0*(v(i,j)+v(i+1,j))*dx)
        df = 0.5d0*(u(i+1,j)-u(i-1,j))*dy+0.5d0*(v(i,j)+v(i+1,j)-v(i,j-1)-v(i+1,j-1))*dx
        ap = aw+ae+as+an+df

        un(i,j) = u(i,j)+dt/dx/dy*(-ap*u(i,j)+aw*u(i-1,j)+ae*u(i+1,j)+as*u(i,j-1)+an*u(i,j+1))-dt*(p(i+1,j)-p(i,j))/dx

        return
        end subroutine upbound_u


!!! compute interior region boundary with 1st-order upwind discrete scheme-->vn
        subroutine upbound_v(N,M,dx,dy,dt,Ra,Pr,u,v,p,T,vn,i,j)
        implicit none
        integer :: N, M, i, j
        real(8) :: dx, dy, dt, Ra, Pr
        real(8) :: u(N,M+1),v(N+1,M),p(N+1,M+1),T(N+1,M+1),vn(N+1,M)
        real(8) :: aw, ae, as, an, df, ap

        aw = Pr+MAX(0.5d0*(u(i-1,j)+u(i-1,j+1))*dy,0.0d0)
        ae = Pr+MAX(0.0d0,-0.5d0*(u(i,j)+u(i,j+1))*dy)
        as = Pr+MAX(0.5d0*(v(i,j-1)+v(i,j))*dx,0.0d0)
        an = Pr+MAX(0.0d0,-0.5d0*(v(i,j)+v(i,j+1))*dx)
        df = 0.5d0*(u(i,j)+u(i,j+1)-u(i-1,j)-u(i-1,j+1))*dy+0.5d0*(v(i,j+1)-v(i,j-1))*dx
        ap = aw+ae+as+an+df

        vn(i,j) = v(i,j)+dt/dx/dy*(-ap*v(i,j)+aw*v(i-1,j)+ae*v(i+1,j) +as*v(i,j-1)+an*v(i,j+1))-dt*(p(i,j+1)-p(i,j))/dy &
                 + dt*Ra*Pr*(T(i,j)+T(i,j+1))/2.0d0

        return
        end subroutine upbound_v


!!! Solve Continuity Equation
        subroutine calpn(N,M,dx,dy,dt,c2,p,un,vn,pn)
        implicit none
        integer :: N, M, i, j
        real(8) :: p(N+1,M+1), un(N,M+1), vn(N+1,M), pn(N+1,M+1)
        real(8) :: dx, dy, dt, c2

        do i=2,N
            do j=2,M
                pn(i,j) = p(i,j)-dt*c2*(  ( un(i,j)-un(i-1,j) )/dx + ( vn(i,j)-vn(i,j-1) ) /dy  )
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


        subroutine caltn(N,M,dx,dy,dt,T,un,vn,Tn)

        implicit none
        integer :: i, j, N, M
        real(8) :: dx, dy, dt, dTx2, dTy2, dTx1, dTy1
        real(8) :: T(N+1,M+1), Tn(N+1,M+1), un(N,M+1), vn(N+1,M)

       ! Interior points using FTCS Sheme(Conservation form)
        do i=3,N-1
            do j=3,M-1
                dTx2 = (T(i+1,j)-2*T(i,j)+T(i-1,j))/dx/dx
                dTy2 = (T(i,j+1)-2*T(i,j)+T(i,j-1))/dy/dy
                dTx1 = ( (un(i+1,j)+un(i,j))/2.0d0*T(i+1,j)-(un(i-1,j)+un(i-2,j))/2.0d0*T(i-1,j) ) / 2.0d0/dx
                dTy1 = ( (vn(i,j+1)+vn(i,j))/2.0d0*T(i,j+1)-(vn(i,j-1)+vn(i,j-2))/2.0d0*T(i,j-1) ) / 2.0d0/dy
                Tn(i,j) = T(i,j)+dt*(dTx2+dTy2-dTx1-dTy1)
            enddo
        enddo

       !Left and right side boundary(Dirichlet B.C.)
        do j=1,M+1
            Tn(1,j) = 4.0d0/3.0d0
            Tn(2,j) = 2.0d0/3.0d0
            Tn(N,j) = 0.0d0
            Tn(N+1,j) = 0.0d0
        enddo

       !Top and bottom side boundary(Neumann B.C.)
        do i=3,N-1
            Tn(i,2) = (4.0d0*Tn(i,3)-Tn(i,4))/3.0d0
            Tn(i,1) = Tn(i,2)
            Tn(i,M) = (4.0d0*Tn(i,M-1)-Tn(i,M-2))/3.0d0
            Tn(i,M+1) = Tn(i,M)
        enddo

        return
        end subroutine caltn


!!! check convergence
        subroutine convergence(N,M,dt,c2,error,u,v,p,T,un,vn,pn,Tn,itc)
        implicit none
        integer :: N, M, i, j, itc
        real(8) :: dt, c2, error, temp
        real(8) :: u(N,M+1), v(N+1,M), p(N+1,M+1), T(N+1,M+1), un(N,M+1), vn(N+1,M), pn(N+1,M+1), Tn(N+1,M+1)
        real(8) :: erru, errv, errp, errT

        itc = itc+1
        erru = 0.0d0
        errv = 0.0d0
        errp = 0.0d0

        do i=1,N
            do j=1,M+1
                temp = ABS(un(i,j)-u(i,j))
                if(temp.GT.erru) erru = temp
                u(i,j) = un(i,j)
            enddo
        enddo

        do i=1,N+1
            do j=1,M
                temp = ABS(vn(i,j)-v(i,j))
                if(temp.GT.errv) errv = temp
                v(i,j) = vn(i,j)
            enddo
        enddo

        do i=1,N+1
            do j=1,M+1
                temp = ABS(pn(i,j)-p(i,j))/c2
                if(temp.GT.errp) errp = temp
                p(i,j) = pn(i,j)
            enddo
        enddo

        do i=1,N+1
            do j=1,M+1
                temp = ABS(Tn(i,j)-T(i,j))
                if(temp.GT.errT) errT = temp
                T(i,j) = Tn(i,j)
            enddo
        enddo

        error = MAX(errp,MAX(erru,errv))
        if(itc.EQ.1) error = 100.0d0
        if(itc.EQ.2) error = 100.0d0

        if(MOD(itc,500).EQ.0) write(*,*) 'itc=',itc,'    |    error=',error

        return
        end subroutine convergence


!!! compute velocity components u, v and pressure p
        subroutine caluvpt(N,M,u,v,p,T,uc,vc,pc,Tc)
        implicit none
        integer :: N, M, i, j
        real(8) :: u(N,M+1), v(N+1,M), p(N+1,M+1), T(N+1,M+1), uc(N,M), vc(N,M), pc(N,M), Tc(N,M)

        do i=1,N
            do j=1,M
                uc(i,j) = 0.5d0*(u(i,j)+u(i,j+1))
                vc(i,j) = 0.5d0*(v(i,j)+v(i+1,j))
                pc(i,j) = 0.25d0*(p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))
                Tc(i,j) = 0.25d0*(T(i,j)+T(i+1,j)+T(i,j+1)+T(i+1,j+1))
            enddo
        enddo

        return
        end subroutine caluvpt


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
            psi(i,M-1) = 0.25d0*psi(i,M-2)
        enddo


        return
        end subroutine calpsi


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
        subroutine output(N,M,X,Y,uc,vc,psi,pc,Tc,k)
        implicit none
        integer :: N, M, i, j, k
        real(8) :: X(N), Y(M), uc(N,M), vc(N,M), psi(N,M), pc(N,M), Tc(N,M)
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
                write(02,100) X(i), Y(j), uc(i,j), vc(i,j), psi(i,j), pc(i,j), Tc(i,j)
            enddo
        enddo

100     format(2x,10(e12.6,'      '))
101     format('Title="Buoyancy Driven Cavity Flow(Artificial Compressibility Methods)"')
102     format('Variables=x,y,u,v,psi,p,T')
103     format('zone',1x,'i=',1x,i5,2x,'j=',1x,i5,1x,'f=point')

        close(02)
        write(*,*) 'Data export to',filename ,'file!'

        return
        end subroutine output


