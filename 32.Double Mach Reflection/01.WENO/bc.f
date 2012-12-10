      
      subroutine bc(rk,uc)
      include "common.h"
      include "para.h"
c--------------------------------------------------------------------------
c     /* set the boundary condition */    
c     for a fifth order weno scheme. we extend 3 ghost cells on each side.
c     + + + -------------------------+ + +
      
c      1. transform the conserved variables to primitive variables.
c      2. set the ghost cell values according to the inner cells.
c      3. re-transform the primitive variables to the conserved variables.
c--------------------------------------------------------------------------
      dimension uc(-3:nxm,-3:nym,4,0:4)
      integer rk
      dimension q(-md:nxm,-md:nym,nq)
      integer m
      parameter(pi=3.141592654)
c
c   1. transform the conserved variables in computational domain to the 
c      primitive variables.
c
      q=0.

      do i=0,nx
      do j=0,ny
      rho=uc(i,j,1,rk)
      xmt=uc(i,j,2,rk)
      ymt=uc(i,j,3,rk)
      eng=uc(i,j,4,rk)
      
      vx=xmt/rho
      vy=ymt/rho

      pre=gm1*(eng-0.5*rho*(vx**2+vy**2))
      q(i,j,1)=rho
      q(i,j,2)=vx
      q(i,j,3)=vy
      q(i,j,4)=pre
      enddo
      enddo

c-------------start bc -----------------------------
c     2. set the south boundary conditon
      do i=0,nx
*    /* south boundary condition.*/
      xxp=1./6.
      xs=real(i)*dx+0.5*dx

      if(xs.le.xxp) then
      do k=-3,-1
      q(i,k,1)=7.996
      q(i,k,2)=8.25*sqrt(3.)/2.
      q(i,k,3)=-8.25*0.5
      q(i,k,4)=116.50

      enddo

      else
      do k=-3,-1
       m=abs(k)-1
       q(i,k,1)= q(i,m,1)
       q(i,k,2)= q(i,m,2)
       q(i,k,3)=-q(i,m,3)
       q(i,k,4)= q(i,m,4)
       enddo
      endif

c       
c     /* north boundary condition */
      xp=1./6.+(1.+20.*time)/sqrt(3.)
      xs=real(i)*dx+0.5*dx
c
      if(xs .le. xp ) then
      do m1=ny+1,ny+3      
      q(i,m1,1)=7.996
      q(i,m1,2)=8.25*sqrt(3.)/2.
      q(i,m1,3)=-8.25*0.5
      q(i,m1,4)=116.50
      enddo
      else   
      do m1=ny+1,ny+3
      q(i,m1,1)=1.4
      q(i,m1,2)=0.0
      q(i,m1,3)=0.0
      q(i,m1,4)=1.0      
      enddo
      endif
      enddo
c

      do j=0,ny
c     /* east boundary condition */
c       " outflow " 
c       " subsonic outflow bc "
      do i=nx+1,nx+3      
      q(i,j,1)=1.4
      q(i,j,2)=2.*q(i-1,j,2)-q(i-2,j,2)
      q(i,j,3)=2.*q(i-1,j,3)-q(i-2,j,3)
      q(i,j,4)=2.*q(i-1,j,4)-q(i-2,j,4)
      enddo
c
c     /* west boundary condition */
c      " inflow "            
c      " supersonic inflow boundary condition "
c
      do i=-1,-3,-1
      q(i,j,1)=7.996
      q(i,j,2)=8.25*cos(1./6.*pi)
      q(i,j,3)=-8.25*sin(1./6.*pi)
      q(i,j,4)=116.50
      enddo
      enddo

c-------------------end  bc --------------------

c***************************************
c    tranform u to uc
      do i=-md,nxm
      do j=-md,nym
      rho=q(i,j,1)
      vx =q(i,j,2)
      vy =q(i,j,3)
      pre=q(i,j,4)

      uc(i,j,1,rk)=rho
      uc(i,j,2,rk)=rho*vx
      uc(i,j,3,rk)=rho*vy
      uc(i,j,4,rk)=pre/gm1+0.5*rho*(vx**2+vy**2)
      enddo
      enddo


     
      return      
      end
