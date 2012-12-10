       
      subroutine rkt(rk,uc)
      include "common.h"
      include "para.h"
c-----------------------------------------------------------
c     runge-kutta in time.
c     mt=3   ---> 3rd r-k (tvd)
c     mt=4   ---> 4rd r-k (tvd)
c-----------------------------------------------------------
      dimension uc(-md:nxm,-md:nym,nq,0:4)
      integer rk

c   3rd TVD R-K Scheme
      if(mt.eq.3) then
      if(rk.eq.0) then
         do m = 1, nq
         do j = 0, ny
         do i = 0, nx
             uc(i,j,m,1)=uc(i,j,m,0)+delt*rhs(i,j,m)
         enddo
         enddo
         enddo

      else if(rk.eq.1) then
          do m = 1, nq
          do j = 0, ny
          do i = 0, nx
             uc(i,j,m,2)=0.75*uc(i,j,m,0) 
     &                 +0.25*(uc(i,j,m,1)+delt*rhs(i,j,m))
          enddo
          enddo
          enddo
      else
          do m = 1, nq
          do j = 0, ny
          do i = 0, nx
            uc(i,j,m,0)=(uc(i,j,m,0) 
     &                  +2*(uc(i,j,m,2)+delt*rhs(i,j,m)))/3.
          enddo
          enddo
          enddo
      endif
      endif

c   4th TVD R-K Scheme
      if(mt.eq.4) then
      if(rk.le.1) then
          do m = 1, nq
          do j = 0, ny
          do i = 0, nx
              uc(i,j,m,rk+1)=uc(i,j,m,0)+0.5*delt*rhs(i,j,m)
          enddo
          enddo
          enddo
      else if(rk.eq.2) then
          do m = 1, nq
          do j = 0, ny
          do i = 0, nx
              uc(i,j,m,3) =uc(i,j,m,0)+delt*rhs(i,j,m)
          enddo
          enddo
          enddo
      else
          do m = 1, nq
          do j = 0, ny
          do i = 0, nx
              uc(i,j,m,0)=(-uc(i,j,m,0)+uc(i,j,m,1) 
     &                    +2*uc(i,j,m,2)+uc(i,j,m,3)
     &                    +0.5*delt*rhs(i,j,m))/3.
          enddo
          enddo
          enddo
      endif
      endif


      return
      end subroutine rkt
