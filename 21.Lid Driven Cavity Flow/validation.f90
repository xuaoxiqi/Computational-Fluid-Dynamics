!!!    This program validate numerical results with reference.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>


        program main
        implicit none
        integer, parameter :: N=129, M=129
        integer :: i, j
        integer :: mid_x, mid_y
        real(8) :: dx, dy
        real(8) :: psi_mid, u_max, v_max, u_max_loc , v_max_loc
        real(8) :: X(N), Y(M), u(N,M), v(N,M), psi(N,M), p(N,M)

        open(unit=01,file='./0041cavity.dat',status='old')
        read(01,*)
        read(01,*)
        read(01,*)
        do j=1,M
            do i= 1,N
                read(01,*) X(i), Y(j), u(i,j), v(i,j), psi(i,j), p(i,j)
            enddo
        enddo
        close(01)

        call validation(N,M,X,Y,u,v,psi,p)

        stop
        end program main


!!! validate results with reference
        subroutine validation(N,M,X,Y,u,v,psi,p)
        implicit none
        integer :: N, M, i, j
        integer :: mid_x, mid_y, temp, temp_max, temp_min
        real(8) :: psi_mid, u_max, v_max, u_max_loc , v_max_loc
        real(8) :: u(N,M), v(N,M), psi(N,M), p(N,M) ,X(N), Y(M)

        mid_x = (N-1)/2+1
        mid_y = (M-1)/2+1
        psi_mid = psi(mid_x,mid_y)

        u_max = 0.0d0
        v_max = 0.0d0


        open(unit=02,file='./u-y.dat',status='unknown')
        open(unit=03,file='./x-v.dat',status='unknown')
        write(02,101)
        write(02,202)
        write(03,302)
        write(02,203) M
        write(03,303) N

        do j=1,M
                write(02,100)  u(mid_x,j), Y(j)
        enddo

        do i=1,N
                write(03,100)  X(i), v(i,mid_y)
        enddo

100     format(2x,10(e12.6,'      '))
101     format('Title="Lid Driven Cavity Flow"')
202     format('Variables=U,Y')
302     format('Variables=X,V')
203     format('zone',1x,'i=',1x,i5,2x,'f=point')
303     format('zone',1x,'i=',1x,i5,2x,'f=point')

        close(02)
        close(03)

        return
        end subroutine validation

