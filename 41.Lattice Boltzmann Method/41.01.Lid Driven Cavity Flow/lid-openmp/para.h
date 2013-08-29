
    integer, parameter :: nx = 129, ny = 129
    real(8), parameter :: dx = 1.0d0/(nx-1), dy = 1.0d0/(ny-1)
    real(8), parameter :: dt = dx

    real(8), parameter :: Re = 100.0d0
    real(8), parameter :: cs2 = 1.0d0/3.0d0
    real(8), parameter :: U_ref = 0.1d0
    real(8), parameter :: tau = 3.0d0*U_ref/Re/dt+0.5d0

    integer, parameter :: itc_max = 500
    real(8), parameter :: eps = 1e-5

