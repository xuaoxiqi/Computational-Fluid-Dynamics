
    integer, parameter :: nx = 129, ny = 129
    real(8), parameter :: dx = 1.0d0/(nx-1), dy = 1.0d0/(ny-1)
    real(8), parameter :: dt = dx

    real(8), parameter :: Ra = 1000.0d0
    real(8), parameter :: Pr = 0.7d0
    real(8), parameter :: Ma = 0.1d0
    real(8), parameter :: cs2 = 1.0d0/3.0d0
    real(8), parameter :: tau_f = 0.5d0+Ma*SQRT(3.0d0*Pr)/dt/SQRT(Ra)
    real(8), parameter :: tau_T = 0.5d0+Ma*2.0d0/dt/SQRT(3.0d0*Ra*Pr)

    integer, parameter :: itc_max = 5*1e5
    real(8), parameter :: eps = 1e-5

