

c     inviscid flux            
      subroutine fx(rk,uc)
c
      include "common.h"
      include "para.h"
      integer nstrt,nend
      dimension w(-md:nxm),vx(-md:nxm),vy(-md:nxm),
     &      h(-md:nxm)
      dimension uc(-md:nxm,-md:nym,4,0:4)
c      real*8  evr(-1:nx,4,4),evl(-1:nx,4,4)
      integer rk,i,j
c 
c
c      write(*,*) 'nxm=',nxm
      em=1.0e-15
c     
c-----begin of outer loop in y-dir-----------------
      do j=0,ny

      am(1)=1.e-15
      am(2)=1.e-15
      am(3)=1.e-15
      am(4)=1.e-15

c     /* compute flux and eigenvalues */
      do i=-md,nxm
      den=uc(i,j,1,rk)
      xmt=uc(i,j,2,rk)
      ymt=uc(i,j,3,rk)
      eng=uc(i,j,4,rk)
      t0=1./den
      vex=xmt*t0
      vey=ymt*t0
      pre=gm1*(eng-0.5*den*(vex*vex+vey*vey))
      
c      write(*,*) 'i=',i
c     write(*,*) 'pressure=', pre 
       ar =sqrt(gamma*pre/den)    
c     
c    write(*,*) 'xmt',xmt
      f(i,1)=xmt
       f(i,2)=vex*xmt+pre
       f(i,3)=vex*ymt
       f(i,4)=vex*(pre+eng)

       uu(i,1)=den
       uu(i,2)=xmt
       uu(i,3)=ymt
       uu(i,4)=eng
      
c      write(*,*)'i=',i,'density=',den
c
      w(i)=sqrt(abs(den))
      h(i)=(pre+eng)*t0
      vx(i)=vex
      vy(i)=vey

      am(1)=max(am(1), abs(vex-ar))
      am(2)=max(am(2), abs(vex))
      am(3)=max(am(3), abs(vex))
      am(4)=max(am(4), abs(vex+ar))

      enddo

      am(1)=1.0*am(1)
      am(2)=1.0*am(2)
       am(3)=1.0*am(3)
      am(4)=1.0*am(4)

      em=max( em, max(am(1),am(4)))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  compute the left and right eigenvectors of roe's mean matrix.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      do i=-1,nx

            t0 = w(i) / ( w(i) + w(i+1) )
            t1 = 1. - t0
          vxm = t0 * vx(i) + t1 * vx(i+1)
          vym = t0 * vy(i) + t1 * vy(i+1)
           hm = t0 *  h(i) + t1 *  h(i+1)
           qm = 0.5 * ( vxm * vxm + vym * vym )
           cm = sqrt( gm1 * ( hm - qm ) )
           t0 = vxm * cm
          evr(i,1,1) = 1.0
          evr(i,1,2) = 0.0
          evr(i,1,3) = 1.0
          evr(i,1,4) = 1.0
          evr(i,2,1) = vxm - cm
          evr(i,2,2) = 0.0
          evr(i,2,3) = vxm
          evr(i,2,4) = vxm + cm
          evr(i,3,1) = vym
          evr(i,3,2) = 1.0
          evr(i,3,3) = vym
          evr(i,3,4) = vym
          evr(i,4,1) = hm - t0
          evr(i,4,2) = vym
          evr(i,4,3) = qm
          evr(i,4,4) = hm + t0
            rcm = 1. / cm
            b1 = gm1 * rcm * rcm
            b2 = qm * b1
            t0 = vxm * rcm
            t1 = b1 * vxm
            t2 = 0.5 * b1
            t3 = b1 * vym
          evl(i,1,1) = 0.5 * ( b2 + t0 )
          evl(i,1,2) = -0.5 * ( t1 + rcm )
          evl(i,1,3) = -0.5 * t3
          evl(i,1,4) = t2
          evl(i,2,1) = - vym
          evl(i,2,2) = 0.0
          evl(i,2,3) = 1.0
          evl(i,2,4) = 0.0
          evl(i,3,1) = 1. - b2
          evl(i,3,2) = t1
          evl(i,3,3) = t3
          evl(i,3,4) = -b1
          evl(i,4,1) =  0.5 * ( b2 - t0 )
          evl(i,4,2) = -0.5 * ( t1 - rcm )
          evl(i,4,3) = -0.5 *   t3
          evl(i,4,4) = t2
      enddo

c     use weno scheme to compute f(i+1/2) and f(i-1/2) 

      call weno_lf(nx)
c
 
      do m=1,nq
       do i=0,nx
       rhs(i,j,m)=(fh(i-1,m)-fh(i,m))*cdx
       enddo
      enddo
 
c      write(*,*) 'j=',j    
      enddo
c-----end loop of y-dir----------------------------
      return
      end


c
      subroutine gy(rk,uc)
      include "common.h"
      include "para.h"
c---------
c      flux g.
c---------      
      dimension w(-md:nym),vx(-md:nym),vy(-md:nym),
     & h(-md:nym)
           dimension uc(-md:nxm,-md:nym,4,0:4)
c      real*8  evr(-1:ny,4,4),evl(-1:ny,4,4)
      integer rk,i,j
 
c--------start loop of x-dir----------------      
      do 200 i=0,nx
      
      am(1)=1.e-15
      am(2)=1.e-15
      am(3)=1.e-15
      am(4)=1.e-15

      do 201 j=-md,nym
      den=uc(i,j,1,rk)
      xmt=uc(i,j,2,rk)
      ymt=uc(i,j,3,rk)
      eng=uc(i,j,4,rk)
      
c      if(den .eq. 0) then
c      write(*,*) 'i=',i
c      endif

       t0=1./den
c
      vex=xmt*t0
      vey=ymt*t0
      pre=gm1*(eng-0.5*den*(vex*vex+vey*vey))
      pre=abs(pre)
      ar=sqrt(gamma*pre*t0)
      
      uu(j,1)=den
      uu(j,2)=xmt
      uu(j,3)=ymt
      uu(j,4)=eng

      f(j,1)=ymt
      f(j,2)=vey*xmt
      f(j,3)=vey*ymt+pre
      f(j,4)=vey*(pre+eng)
      
      den=abs(den)
      w(j)=sqrt(den)
      h(j)=(pre+eng)*t0
      vx(j)=vex
      vy(j)=vey
      am(1)=max(am(1),abs(vey-ar))
      am(2)=max(am(2),abs(vey))
      am(3)=max(am(3),abs(vey))
      am(4)=max(am(4),abs(vey+ar))
  201 continue
c
c     a small trick 
c      am(1)=am(1)*1.1
c       am(2)=am(2)*1.1
c        am(3)=am(3)*1.1
c       am(4)=am(4)*1.1

c      em=max(em,max(am(1),am(4)))
c
c     compute the left and right eigenvectors of roe's mean matrix.
      do 202 j=-1,ny

          t0 = w(j) / ( w(j) + w(j+1) )
          t1 = 1. - t0
          vxm = t0 * vx(j) + t1 * vx(j+1)
          vym = t0 * vy(j) + t1 * vy(j+1)
           hm = t0 *  h(j) + t1 *  h(j+1)
           qm = 0.5 * ( vxm * vxm + vym * vym )
           cm = sqrt( gm1 * ( hm - qm ) )
           t0 = vym * cm
          evr(j,1,1) = 1.0
          evr(j,1,2) = 0.0
          evr(j,1,3) = 1.0
          evr(j,1,4) = 1.0
          evr(j,2,1) = vxm
          evr(j,2,2) = 1.0
          evr(j,2,3) = vxm
          evr(j,2,4) = vxm
          evr(j,3,1) = vym - cm
          evr(j,3,2) = 0.0
          evr(j,3,3) = vym
          evr(j,3,4) = vym + cm
          evr(j,4,1) = hm - t0
          evr(j,4,2) = vxm
          evr(j,4,3) = qm
          evr(j,4,4) = hm + t0
            rcm = 1. / cm
            b1 = gm1 * rcm * rcm
            b2 = qm * b1
            t0 = vym * rcm
            t1 = b1 * vym
            t2 = 0.5 * b1
            t3 = b1 * vxm
          evl(j,1,1) = 0.5 * ( b2 + t0 )
          evl(j,1,2) = -0.5 * t3
          evl(j,1,3) = -0.5 * ( t1 + rcm )
          evl(j,1,4) = t2
          evl(j,2,1) = - vxm
          evl(j,2,2) = 1.0
          evl(j,2,3) = 0.0
          evl(j,2,4) = 0.0
          evl(j,3,1) = 1. - b2
          evl(j,3,2) = t3
          evl(j,3,3) = t1
          evl(j,3,4) = -b1
          evl(j,4,1) =  0.5 * ( b2 - t0 )
          evl(j,4,2) = -0.5 *   t3
          evl(j,4,3) = -0.5 * ( t1 - rcm )
          evl(j,4,4) = t2
  202 continue
      
      call weno_lf(ny)
      
      do m=1,nq
      do j=0,ny
      rhs(i,j,m)=rhs(i,j,m)+(fh(j-1,m)-fh(j,m))*cdy
      enddo
      enddo     
  200 continue
c---------end loop ----------------
      return
      end
