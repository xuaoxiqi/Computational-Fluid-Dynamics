	integer, parameter :: nx = 129, ny = 129
	real(8), parameter :: Re = 1000.0d0
	real(8), parameter :: dt = 0.0001d0
	real(8), parameter :: dx = 1.0d0/(nx-1), dy = 1.0d0/(ny-1)
	real(8), parameter :: eps = 1e-6
	integer, parameter :: itc_max = 2000
