c-------------------------------------------------
c      output the final solutions.
c
c-------------------------------------------------
ccccccccccccccccccccccccc

	subroutine output(uc,r,rk)
      include "para.h"
	include "common.h"
	dimension uc(-md:nxm,-md:nym,4,0:3)
      dimension r(0:nx,0:ny,2)
      character (len=100):: filename
      integer i,j
      integer rk
      
      write(filename,*)nstep
	filename = adjustl(filename)
      open(02,file='output_'//trim(filename)//'.plt',status='unknown')
      	write(02,201)
      write(02,202)
      write(02,203) 1+nx, 1+ny
 	do j=0,ny
 	    do i=0,nx
              rho=uc(i,j,1,rk)
	        vx=uc(i,j,2,rk)/rho
	        vy=uc(i,j,3,rk)/rho
	        p=gm1*(uc(i,j,4,rk)-0.5*rho*(vx**2+vy**2))

 	        write(02,200) r(i,j,1),r(i,j,2),rho,vx,vy,p
 	    enddo
  	enddo

200   format(2x,10(e12.6,'      '))
201   format('Title="Double Mach Reflection Problem"')
202   format('Variables=x,y,rho,u,v,p')
203   format('zone',1x,'i=',1x,i5,2x,'j=',1x,i5,1x,'f=point')

      close(02)

      return
	end subroutine output