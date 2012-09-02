!!!    This program validate lid driven cavity numerical results with reference.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>


        program main
        implicit none
        integer, parameter :: N=81, M=81
        integer :: i, j
        integer :: mid_x, mid_y, temp, temp_max, temp_min
        real(8) :: dx, dy
        real(8) :: psi_mid, u_max, v_max, u_max_loc , v_max_loc
        real(8) :: X(N), Y(M), u(N,M), v(N,M), psi(N,M)

        open(unit=01,file='./ra6.dat',status='old')
        read(01,*)
        read(01,*)
        read(01,*)
        do j=1,M
            do i= 1,N
                read(01,*) X(i), Y(j), u(i,j), v(i,j), psi(i,j)
            enddo
        enddo
        close(01)

        dx = 1.0d0/(N-1)
        dy = 1.0d0/(M-1)

        open(unit=03,file='results.txt',status='unknown')

        call validation(N,M,dx,dy,u,v,psi)

        close(03)

        stop
        end program main


!!! validate results with reference
        subroutine validation(N,M,dx,dy,u,v,psi)
        implicit none
        integer :: N, M, i, j
        real(8) :: dx, dy
        real(8) :: psi_max
        real(8) :: u(N,M), v(N,M), psi(N,M)

        psi_max = 0.0d0

        do i=1,N
            do j=1,M
                if(psi(i,j).GT.psi_max) then
                    psi_max = psi(i,j)
                endif
            enddo
        enddo


        write(03,*)
        write(03,*) 'psi_max =',psi_max
        write(03,*)

        return
        end subroutine validation

