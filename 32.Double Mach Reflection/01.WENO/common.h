      
c
      common /a/ me,mt 
	integer me,mt 

c     boundary type.
	common/b/  bce,bcs,bcw,bcn
	integer bce,bcs,bcw,bcn
c
	common/domain/ nxm,nym
	integer nxm,nym

c       
c	common/len/ lx,ly
c	real*8 lx,ly
      
	common/del/ dx,dy,cdx,cdy
	real*8 dx,dy,cdx,cdy

      common/para/cfl,tend,fsmach,pr,
     &  delt,em,am(5)

c     subsonic inflow/outflow bc conditions.  initialize to zero.  if input
c     does not reset to non-zero value, then set based on input freestream.
c     input t_total as kelvins or degrees r,
c           p_back, p_total as p/p_inf
c     and non-dimensionalize correctly at end of routine.
c     
      common/bcio/ pbakbc,ptotbc,ttotbc
      real*8  pbakbc,ptotbc,ttotbc
	
	common/free/ rhoinf,ainf 
      
	common/rhs/ rhs(0:500,0:500,4)
c
	common/ff/ f(-3:500,4),uu(-3:500,4),fh(-3:500,4)
	real*8 f,uu

      common/time/ time,nstep
	real*8 time
      integer nstep
c

	common/eigv/ evr(-1:500,4,4),evl(-1:500,4,4)
	real*8 evr,evl


c
