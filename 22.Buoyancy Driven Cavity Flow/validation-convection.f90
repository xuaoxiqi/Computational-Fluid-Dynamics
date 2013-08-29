!!!    This program validate numerical results with benchmark reference.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>


    program main
    implicit none
    integer, parameter :: nx=81, ny=81
    integer :: i, j
    real(8) :: dx, dy
    real(8) :: X(nx), Y(ny), u(nx,ny), v(nx,ny), T(nx,ny)

    X = 0.0d0
    Y = 0.0d0
    T = 0.0d0
    u = 0.0d0
    v = 0.0d0

    dx = 1.0d0/(nx-1)
    dy = 1.0d0/(ny-1)

    open(unit=01,file='./ra1e3.plt',status='old')
    read(01,*)
    read(01,*)
    read(01,*)
    do j=1,ny
        do i= 1,nx
            read(01,*) X(i), Y(j), T(i,j), u(i,j), v(i,j)
        enddo
    enddo
    close(01)

    call validation(nx,ny,dx,dy,u,v,T)

    stop
    end program main


!!! validate results with reference
    subroutine validation(nx,ny,dx,dy,u,v,T)
    implicit none
    integer :: nx, ny, i, j
    integer :: mid_x, mid_y, temp, temp_max, temp_min
    real(8) :: dx, dy
    real(8) :: u_max, v_max, u_max_loc , v_max_loc
    real(8) :: Nu_0, Nu_max, Nu_min, Nu_max_loc, Nu_min_loc
    real(8) :: u(nx,ny), v(nx,ny), T(nx,ny)
    real(8) :: Nu(ny)

    mid_x = (nx+1)/2
    mid_y = (ny+1)/2

    u_max = 0.0d0
    v_max = 0.0d0

    temp = 0
    do j=1,ny
        if(u(mid_x,j).GT.u_max) then
            u_max = u(mid_x,j)
            temp = j
        endif
    enddo
    u_max_loc = (temp-1)*dy

    temp = 0
    do i=1,nx
        if(v(i,mid_y).GT.v_max) then
            v_max = v(i,mid_y)
            temp = i
        endif
    enddo
    v_max_loc = (temp-1)*dx

    temp_max = 0
    temp_min = 0
    Nu_max = 0.0d0
    Nu_min = 100.0d0
    Nu_0 = 0.0d0
    do j=1,ny
        !!Nu(j) = -(-11.0d0*T(1,j)+18.0d0*T(2,j)-9.0d0*T(3,j)+2.0d0*T(4,j))/6.0d0/dy
        Nu(j) = -(-3.0d0*T(1,j)+4.0d0*T(2,j)-T(3,j))/2.0d0/dy
        !!Nu(j) = -(T(2,j)-T(1,j))/dx
        Nu_0 = Nu_0+Nu(j)
        if(Nu(j).GT.Nu_max) then
            Nu_max = Nu(j)
            temp_max = j
        elseif(Nu(j).LT.Nu_min) then
            Nu_min = Nu(j)
            temp_min = j
        endif
    enddo
    Nu_0 = Nu_0/ny
    Nu_max_loc = (temp_max-1)*dy
    Nu_min_loc = (temp_min-1)*dy

    open(unit=03,file="validation.txt",status="unknown")
    write(03,*) "1: Maximum velocity on mid-plane."
    write(03,*) "U_max=",u_max,"; at x=0.5, y= ",u_max_loc
    write(03,*) "V_max=",v_max,"; at y=0.5, x= ",v_max_loc
    write(03,*) " "
    write(03,*) "2: Nusselt number on x = 0 plane."
    write(03,*) "Nu_0 = ", Nu_0
    write(03,*) "Nu_max=",Nu_max,"; at y= ",Nu_max_loc
    write(03,*) "Nu_min=",Nu_min,"; at y = ",Nu_min_loc
    close(03)

    return
    end subroutine validation

