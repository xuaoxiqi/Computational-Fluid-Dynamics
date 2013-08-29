    common /flow/ X, Y, u, v, rho, p, psi, un
    real(8) :: X(nx), Y(ny)
    real(8) :: u(nx,ny), v(nx,ny), rho(nx,ny), p(nx,ny), psi(nx,ny)
    real(8) :: un(0:8)

    common /lbm/ omega, f
    real(8) :: omega(0:8)
    real(8) :: f(0:8,nx,ny)

    common /para/ itc
    integer :: itc

    real(8) :: ex(0:8), ey(0:8)
    data ex/0.0d0,1.0d0,0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0/
    data ey/0.0d0,0.0d0,1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0/

