    common /flow/ X, Y, u, v, p, psi
    real(8) :: X(nx), Y(ny)
    real(8) :: u(nx,ny+1), v(nx+1,ny), p(nx+1,ny+1), psi(nx,ny)

    common /para/ itc
    integer :: itc

    common /guess/ uGuess, vGuess, pGuess
    real(8) :: uGuess(nx,ny+1), vGuess(nx+1,ny), pGuess(nx+1,ny+1)
