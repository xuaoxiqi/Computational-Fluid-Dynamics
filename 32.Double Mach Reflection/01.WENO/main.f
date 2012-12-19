cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccW
c        two dimensional weno schemes.
c        uniform grid.
c        apply multiple boundary conditions
c        this program is used to solve n-s eqaution.
c        u_t= rhs= -f(u)_x-f(u)_y+fv(v)_x+gv(u)_y.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      program main
      include "para.h"
      include "common.h"
      dimension r(0:500,0:200,2)
      dimension  q(-3:503,-3:153,4)
     &    , qc(-3:503,-3:153,4,0:4)
      integer rk
    
c     /* read parameters */
      call setup

c     /* initialization */
      call init(r,q,qc)    
      call wrinit(r,q)
      
      rk=0
      delt=cfl*max(dx,dy)
      nstep=0

      write(*,*) 'initial delta t=',delt

      do    
          call bc(rk,qc)
          call fx(rk,qc)
          call gy(rk,qc)
          call rkt(rk,qc)

          rk=mod(rk+1,mt)

          if(rk.eq.0) then
              time=time+delt
              write(*,*) 'time=',time
          endif
      
          if(time.gt.tend) exit
          nstep=nstep+1
          write(*,*) 'steps=',nstep

          if( mod(nstep,100).eq. 0) then
                    call output(qc,r,rk)
          endif

      enddo

      rk = 0
      call output(qc,r,rk)
      
      write(*,*) 'End!'
      
      stop      
      end program main
