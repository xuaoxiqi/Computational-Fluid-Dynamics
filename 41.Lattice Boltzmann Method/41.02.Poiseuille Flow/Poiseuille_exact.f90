

        program main
        implicit none

        integer, parameter :: N=17,M=17
        integer :: mid_x, mid_y
        integer :: i, j
        real(8) :: u_max
        real(8) :: nu
        real(8) :: dx, dy
        real(8) :: X(N), Y(M)

        nu = 0.01d0

        u_max = 0.001/8.0d0/nu

        write(*,*) "u_max",u_max

        dx = 1.0d0/N
        dy = 1.0d0/M
        do i=1,N
            X(i) = (i-1)*dx
        enddo
        do j=1,M
            Y(j) = (j-1)*dy
        enddo

        open(unit=02,file='exact.dat',status='unknown')
        write(02,101)
        write(02,202)
        write(02,203) M-1


        do j=1,M-1
                write(02,100)  4.0d0*u_max/M/M*j*(M-j), Y(j+1)
        enddo


100     format(2x,10(e12.6,'      '))
101     format('Title="Lid Driven Cavity Flow"')
202     format('Variables=U,Y')
203     format('zone',1x,'i=',1x,i5,2x,'f=point')

        close(02)

        stop
        end program main
