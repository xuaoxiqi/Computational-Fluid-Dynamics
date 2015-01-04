
        ! predictor
        do i=1,nx-1
            uPre(i) = u(i)-constC*dt/dx*(u(i+1)-u(i))+nu*dt/dx/dx*(u(i+1)-2.0d0*u(i)+u(i-1))
        enddo
        uPre(0) = 100.0d0
        uPre(nx) = 0.0d0

        ! corrector
        do i=1,nx-1
            u(i) = 0.5d0*( u(i)+uPre(i)-constC*dt/dx*(uPre(i)-uPre(i-1))+nu*dt/dx/dx*(uPre(i+1)-2.0d0*uPre(i)+uPre(i-1)) )
        enddo

        u(0) = 100.0d0
        u(nx) = 0.0d0
