      
      subroutine weno_lf(ns)
      include "common.h"
      include "para.h"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this routine is used to construct weno scheme. 
c     5-th order weno scheme
c     use local characteristic decomposition.
c         lax-friedriches flux splitting technique.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      integer ns,i,j
      integer nsm
c    f: inviscid flux 
c       f(u) or g(u)
c    u: the conserved variables.
c  evl: left eig-vector
c  evr: right eig-vector
c   fh: flux at half nodes.
*

c     local variables:
      dimension ff(-md:ns+md,nq)
      dimension gg(-md:ns+md,nq,2)
      dimension hh(-md:ns+md,nq,2)
      dimension df(-md:ns+md,nq)
      dimension du(-md:ns+md,nq)
        double precision :: tau5
c
      nsm=ns+md-1

c   compute the difference of variables & flux
      do m=1,nq
      do i=-md,nsm-1
        df(i,m)=f(i+1,m)-f(i,m)
        du(i,m)=uu(i+1,m)-uu(i,m)
      enddo
      enddo

c------------loop over in " m" characteristic field-------------------
      do 100 m=1 ,nq
c
c       *  use lax-friedriches method to split the fluxes.
      do m1=1,nq
      do i=-md,nsm-1
      gg(i,m1,1)=0.5*(df(i,m1)+am(m)*du(i,m1))
      gg(i,m1,2)=gg(i,m1,1)-df(i,m1)
      enddo
      enddo
 
c
c    /*   fifth order */

c    project the positive and negative part part of the fluxes the
c      'm' th charactersitic field.
      
      do m1=1,nq
      k0=m1-3
      k1=3-m1
      do i=-1,ns
          hh(i,m1,1) = evl(i,m,1)*gg(i+k0,1,1) + evl(i,m,2)*gg(i+k0,2,1)
     &               + evl(i,m,3)*gg(i+k0,3,1) + evl(i,m,4)*gg(i+k0,4,1)
          hh(i,m1,2) = evl(i,m,1)*gg(i+k1,1,2) + evl(i,m,2)*gg(i+k1,2,2)
     &               + evl(i,m,3)*gg(i+k1,3,2) + evl(i,m,4)*gg(i+k1,4,2)
      enddo
      enddo

c* compute the weights and approximate the fluxes

      do i = -1, ns
         ff(i,m) = 0.0
      enddo

      do m1 = 1, 2
      do i = -1, ns
          t1 = hh(i,1,m1) - hh(i,2,m1)
          t2 = hh(i,2,m1) - hh(i,3,m1)
          t3 = hh(i,3,m1) - hh(i,4,m1)

c      jiang & shu

          tt1 = 13. * t1**2 + 3. * (   hh(i,1,m1) - 3*hh(i,2,m1) )**2
          tt2 = 13. * t2**2 + 3. * (   hh(i,2,m1) +   hh(i,3,m1) )**2
          tt3 = 13. * t3**2 + 3. * ( 3*hh(i,3,m1) -   hh(i,4,m1) )**2

c          tt1 =  ( eps + tt1 )**2
c          tt2 =  ( eps+ tt2 )**2
c          tt3 =  ( eps + tt3 )**2
 
c---------improved weno------------------
           p=2.d0
c
           tau5=abs(tt1-tt3)
           s1=0.1d0*(1.d0+(tau5/(tt1+eps))**p)
           s2=0.6d0*(1.d0+(tau5/(tt2+eps))**p)
           s3=0.1d0*(1.d0+(tau5/(tt3+eps))**p)        

            t0 = 1. / ( s1 + s2 + s3 )
            s1 = s1 * t0
            s3 = s3 * t0
          ff(i,m) = ff(i,m) + ( s1*(t2-t1) + (0.5*s3-0.25)*(t3-t2) ) /3.
      enddo
      enddo

  100 continue    
c------------ end loop------------------------------------------------

c project the fluxes to the physical space.

      
      do m =  1, nq
      do i = -1, ns
         fh(i,m) =evr(i,m,1) * ff(i,1) + evr(i,m,2) * ff(i,2)
     &           +evr(i,m,3) * ff(i,3) + evr(i,m,4) * ff(i,4)
     &           +(-f(i-1,m)+7*(f(i,m)+f(i+1,m))-f(i+2,m))/12.
      enddo
      enddo

      return
      end
