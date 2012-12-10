      
c    set up for the parameters.
c  
      subroutine setup
      include "common.h"
      include "para.h"

      me=1    !Equation type: 1->Euler;
      mt=4    !Order in time
      nxm=nx+md
      nym=ny+md
      cfl=0.05d0
      tend=0.20
     
      write(*,*) '************************'
      write(*,*) 'Double Mach Reflection Problem'
      write(*,*) 'Order in time',mt
      if( me.eq.1) then
          write(*,*) 'Euler Equation'
      else
          write(*,*) 'Navier Stokes Equation'
      endif
      
      write(*,*) 'cfl=',cfl
      write(*,*) '*************************'
c
c
      return
      end subroutine setup


   
      subroutine init(r,u,uc)
      include "common.h"
      include "para.h"
      dimension r(0:nx,0:ny,2)
      dimension u(-md:nx+md,-md:ny+md,nq)
      dimension uc(-md:nx+md,-md:ny+md,nq,0:4)
      integer i,j

c     grid information.
      xleft=0.0d0
      xright=xleft+lx
      ydown=0.0d0
      yup=ydown+ly

      dx=(xright-xleft)/(nx+1)
      cdx=1.0d0/dx

      dy=(yup-ydown)/(ny+1)
      cdy=1.0d0/dy

c     nodes information
      do j=0,ny
          do i=0,nx
              r(i,j,1)=i*dx+xleft+0.5*dx
              r(i,j,2)=j*dy+ydown+0.5*dy
          enddo
      enddo

      vx= 8.25d0*SQRT(3.0d0)/2.0d0
      vy=-8.25d0*0.5d0
      
      prel=116.5d0
      prer=1.0d0

      rhol=7.9996d0
      rhor=1.4d0
       
      x0=1.0d0/6.0d0

c... initialize the condition.
      do j=0,ny
        do i=0,nx
            x=r(i,j,1)
            y=r(i,j,2)
         
c         xcut=x0+y*(1./3)**(0.5)
            xcut =x0+y*1./sqrt(3.)
            if(x.le.xcut) then
                u(i,j,1)=rhol
                u(i,j,2)=vx
                u(i,j,3)=vy
                u(i,j,4)=prel
            else
                u(i,j,1)=rhor
                u(i,j,2)=0.
                u(i,j,3)=0.
                u(i,j,4)=prer
            endif
        enddo
      enddo
      
      do j=0,ny
        do i=0,nx
            rho=u(i,j,1)
            vx=u(i,j,2)
            vy=u(i,j,3)
            pre=u(i,j,4)

            uc(i,j,1,0)=rho
            uc(i,j,2,0)=rho*vx
            uc(i,j,3,0)=rho*vy
            uc(i,j,4,0)=pre/gm1+rho*(vx*vx+vy*vy)/2.d0
         enddo
      enddo
  
      return
      end subroutine init

c
      subroutine wrinit(r,u)
      include "common.h"
      include "para.h"
      dimension r(0:nx,0:ny,2)
      dimension u(-md:nx+md,-md:ny+md,4)
      integer i,j
      
      open(01,file='init.plt',status='unknown')
      write(01,101)
      write(01,102)
      write(01,103) 1+nx, 1+ny
      do j=0,ny
          do i=0,nx
              write(01,100) r(i,j,1),r(i,j,2),u(i,j,1),u(i,j,2),
     &           u(i,j,3),u(i,j,4)
      enddo
      enddo

100   format(2x,10(e12.6,'      '))
101   format('Title="Double Mach Reflection Problem"')
102   format('Variables=x,y,rho,u,v,p')
103   format('zone',1x,'i=',1x,i5,2x,'j=',1x,i5,1x,'f=point')

      close(01)
      
      return
      end subroutine wrinit
