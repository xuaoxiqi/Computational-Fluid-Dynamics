!!!    This program validate numerical results with reference.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>


        program main
        implicit none
        integer, parameter :: N=81, M=81
        integer :: i, j
        integer :: mid_x, mid_y, temp, temp_max, temp_min
        real(8) :: dx, dy
        real(8) :: psi_mid, u_max, v_max, u_max_loc , v_max_loc, Nu_max, Nu_min, Nu_max_loc, Nu_min_loc
        real(8) :: X(N), Y(M), u(N,M), v(N,M), psi(N,M), T(N,M)
        real(8) :: Nu(N)

        open(unit=01,file='./1104cavity.dat',status='old')
        read(01,*)
        read(01,*)
        read(01,*)
        do j=1,M
            do i= 1,N
                read(01,*) X(i), Y(j), u(i,j), v(i,j), psi(i,j), T(i,j)
            enddo
        enddo
        close(01)

        dx = 1.0d0/(N-1)
        dy = 1.0d0/(M-1)

        open(unit=03,file='results.txt',status='unknown')

        call validation(N,M,dx,dy,u,v,psi,T)

        close(03)

        stop
        end program main


!!! validate results with reference
        subroutine validation(N,M,dx,dy,u,v,psi,T)
        implicit none
        integer :: N, M, i, j
        integer :: mid_x, mid_y, temp, temp_max, temp_min
        real(8) :: dx, dy
        real(8) :: psi_mid, u_max, v_max, u_max_loc , v_max_loc, Nu_max, Nu_min, Nu_max_loc, Nu_min_loc
        real(8) :: u(N,M), v(N,M), psi(N,M), T(N,M)
        real(8) :: Nu(N)

        mid_x = INT(N/2)
        mid_y = INT(M/2)
        psi_mid = psi(mid_x,mid_y)

        u_max = 0.0d0
        v_max = 0.0d0
        Nu_max = 0.0d0
        Nu_min = 100.0d0
        temp = 0
        temp_max = 0
        temp_min = 0

        do j=1,M
            if(u(mid_x,j).GT.u_max) then
                u_max = u(mid_x,j)
                temp = j
            endif
        enddo
        u_max_loc = (temp-1)*dy

        do i=1,N
            if(v(i,mid_y).GT.v_max) then
                v_max = v(i,mid_y)
                temp = i
            endif
        enddo
        v_max_loc = (temp-1)*dx

        do j=1,M
            !!!Nu(j) = -(-11.0d0*T(1,j)+18.0d0*T(2,j)-9.0d0*T(3,j)+2.0d0*T(4,j))/6.0d0/dy
            !!!Nu(j) = -(-3.0d0*T(1,j)+4.0d0*T(2,j)-T(3,j))/2.0d0/dy
            Nu(j) = -(T(2,j)-T(1,j))/dx
            if(Nu(j).GT.Nu_max) then
                Nu_max = Nu(j)
                temp_max = j
            elseif(Nu(i).LT.Nu_min) then
                Nu_min = Nu(j)
                temp_min = j
            endif
        enddo
        Nu_max_loc = (temp_max-1)*dx
        Nu_min_loc = (temp_min-1)*dx

        write(03,*)
        write(03,*) 'psi_mid =',psi_mid
        write(03,*) 'u_max =',u_max,'at y =',u_max_loc
        write(03,*) 'v_max =',v_max,'at x =',v_max_loc
        write(03,*) 'Nu_max =',Nu_max,'at y=',Nu_max_loc
        write(03,*) 'Nu_min =',Nu_min,'at y=',Nu_min_loc
        write(03,*)

        return
        end subroutine validation

