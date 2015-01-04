
        do j=2,ny-1
            do i=2,nx-1
                T(i,j) = T(i,j)+dt*alpha*( (T(i+1,j)-2.0d0*T(i,j)+T(i-1,j))/dx/dx+(T(i,j+1)-2.0d0*T(i,j)+T(i,j-1))/dy/dy )
            enddo
        enddo

        ! Left and right side B.C.
        do j=1,ny
            T(1,j) = 0.0d0
            T(nx,j) = 0.0d0
        enddo
        ! Top and bottom side B.C.
        do i=1,nx
            T(i,1) = (4.0d0*T(i,2)-T(i,3))/3.0d0
            T(i,ny) = dsin(Pi*X(i))
        enddo
