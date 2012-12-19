      PROGRAM main
c     This program deals with 1d Euler equations (du/dt+df(u)/dx=0) with weno scheme.
      implicit none
c     const. for 3rd weno
      real(8)::cr(0:2,-1:2)=(/ 11.d0/6.,-7.d0/6., 1.d0/3.,
     $                          1.d0/3., 5.d0/6.,-1.d0/6.,
     $                         -1.d0/6., 5.d0/6., 1.d0/3.,
     $                          1.d0/3.,-7.d0/6., 11.0d0/6./)
      real(8)::dr(1:2,0:2)=(/ 0.3d0,0.1d0,0.6d0,0.6d0,0.1d0,0.3d0 /)
      common /group1/ cr,dr
c     const. for gas properties
      real(8)::gama=1.4d0,cv=1.0d0
      common /group2/ gama,cv
c     variables for domain and grid
      real(8) xmax,xmin,dx,dt,t
      real(8),allocatable::x(:)
      integer N,itert,itc
c     varibles for flow field
      real(8),allocatable::u(:,:),utmp(:,:),flux(:,:)
c     varibles for controling the B.C.
c     xbc=1 ==> periodic in x-direction
      integer xbc,i,j

      xmax = 1.0d0
      xmin = 0.0d0
      N = 100
      dt = 1e-5
      t = 0.25d0
      itert = t/dt


      allocate(x(-3:N+3))
      call grid(N,xmax,xmin,dx,x)
      allocate(u(3,-3:N+3))
      allocate(utmp(3,-3:N+3))
      call initialize(u,N,x)
      call bondc(xbc)
      allocate(flux(3,-1:N))

      write(*,'(a)')'==========WENO Scheme 1D system=========='

      do itc=1,itert
        call weno(u,N,flux,dx)
        do i=0,N
            utmp(:,i)=u(:,i)-dt/dx*(flux(:,i)-flux(:,i-1))
        enddo
        call weno(utmp,N,flux,dx)
        do i=0,N
            utmp(:,i)=0.75*u(:,i)+0.25*utmp(:,i)-0.25*dt/dx
     $             *(flux(:,i)-flux(:,i-1))
        enddo
        call weno(utmp,N,flux,dx)
        do i=0,N
            u(:,i)=1.0d0/3.0*u(:,i)+2.0d0/3.0*utmp(:,i)-2.0d0/3.0*dt/dx
     $             *(flux(:,i)-flux(:,i-1))
        enddo
      enddo
      call output(u,x,N)

      deallocate(x)
      deallocate(u)
      deallocate(utmp)
      deallocate(flux)
      END PROGRAM main
c     ================================================================

c     ================================================================

      subroutine grid(N,xmax,xmin,dx,x)
      real(8) xmax,xmin,dx
      integer N,i
      real(8) x(-3:N+3)

      dx=(xmax-xmin)/N
      do i=-3,N+3
        x(i)=xmin+i*dx   
      enddo
      end
c     =====================================================

c     =====================================================
      subroutine initialize(u,N,x)
      integer N,i
      real(8) u(3,-3:N+3),x(-3:N+3),pi
      real(8) gama,cv
      common /group2/ gama,cv
      pi=acos(-1.0d0)

      do i=0,N
        if(x(i).LE.0.5) then
            u(1,i) = 1.0d0
            u(2,i) = 0.0d0
            u(3,i) = 2.5d0
        else
            u(1,i) = 0.125d0
            u(2,i) = 0.0d0
            u(3,i) = 0.25d0
        endif
      enddo

      end subroutine initialize
c     =====================================================

c     =====================================================
      subroutine bondc(xbc)
      integer xbc

      xbc=1
      end
c     =====================================================

c     =====================================================
      subroutine weno(u,N,flux,dx)
      integer N,ix,i,j
      real(8) u(3,-3:N+3),flux(3,-1:N)
      real(8) dx,gama,cv
      real(8) cr(0:2,-1:2),dr(1:2,0:2)
      common /group1/ cr,dr
      common /group2/ gama,cv

      real(8) lambda(3),leftv(3,3),rightv(3,3),vflux(3)
      real(8) vr(0:2),belta(0:2),alpha(0:2),omiga(0:2)
      real(8) vstencil(-2:2),fu(3),vsv(-2:2)
      real(8) g(3,-2:3),v(3,-2:3),lambdamax(3)
      real(8) ro,p,vel,c,rol,ror,pl,pr,ul,ur,h,hl,hr

      call ghostp(u,N,dx)
c     -------- caculate the flux at x(i+1/2),i=-1,...,N --------
      do ix=-1,N
c     Roe average
        rol=sqrt(u(1,ix))
        ror=sqrt(u(1,ix+1))
        ro=rol*ror
        ul=u(2,ix)/u(1,ix)
        ur=u(2,ix+1)/u(1,ix+1)
        vel=(rol*ul+ror*ur)/(rol+ror)
        pl=(gama-1.d0)*(u(3,ix)  -0.5*u(1,ix)  *ul*ul)
        pr=(gama-1.d0)*(u(3,ix+1)-0.5*u(1,ix+1)*ur*ur)
        p=(rol*pl+ror*pr)/(rol+ror)
        hl=(u(3,ix)  + pl)/u(1,ix)
        hr=(u(3,ix+1)+ pr)/u(1,ix+1)
        h=(hl*rol+hr*ror)/(rol+ror)
        c=sqrt((gama-1.d0)*(h-0.5*vel*vel))
        if(p<0.)then
            write(*,*) "P < 0 at 1 !"
            stop
        endif
        rightv=1.d0
        rightv(2,1)=vel-c
        rightv(3,1)=vel*vel*0.5-vel*c+c*c/(gama-1.d0)
        rightv(2,2)=vel
        rightv(3,2)=vel*vel*0.5
        rightv(2,3)=vel+c
        rightv(3,3)=vel*vel*0.5+vel*c+c*c/(gama-1.d0)

        leftv(1,1)=0.25*vel*vel*(gama-1.d0)/(c*c)+0.5*vel/c
        leftv(2,1)=1.d0-0.5*(gama-1.d0)*vel*vel/(c*c)
        leftv(3,1)=0.25*vel*vel*(gama-1.d0)/(c*c)-0.5*vel/c
        leftv(1,2)=0.5*(1.d0-gama)*vel/(c*c)-0.5/c
        leftv(2,2)=(gama-1.d0)*vel/(c*c)
        leftv(3,2)=0.5*(1.d0-gama)*vel/(c*c)+0.5/c
        leftv(1,3)=0.5*(gama-1.d0)/(c*c)
        leftv(2,3)=(1.d0-gama)/(c*c)
        leftv(3,3)=0.5*(gama-1.d0)/(c*c)

        lambdamax=0.d0
        do i=-2,3
        fu(1)=u(2,i+ix)
        fu(2)=(gama-1.)*u(3,i+ix)+0.5*(3.-gama)
     $       *u(2,i+ix)*u(2,i+ix)/u(1,i+ix)
        fu(3)=(gama*u(3,i+ix)-0.5*(gama-1.)*u(2,i+ix)
     $       *u(2,i+ix)/u(1,i+ix))*u(2,i+ix)/u(1,i+ix)
        g(:,i)=matmul(leftv,fu)
        v(:,i)=matmul(leftv,u(:,i+ix))

        ro=u(1,ix+i)
        vel=u(2,ix+i)/ro
        p=(gama-1.d0)*(u(3,ix+i)-0.5*ro*vel*vel)
        if(p<0.)then
            write(*,*) "P < 0 at 2 !"
            stop
        endif
        c=sqrt(gama*p/ro)
        if(abs(vel-c)>lambdamax(1)) lambdamax(1)=abs(vel-c)
        if(abs(vel  )>lambdamax(2)) lambdamax(2)=abs(vel  )
        if(abs(vel+c)>lambdamax(3)) lambdamax(3)=abs(vel+c)
        enddo

        do i=1,3     
c     -------- caculate the g(i+1/2)+ ------------
            vstencil=g(i,-2:2)
            vsv=v(i,-2:2)
            
            vstencil=0.5*(vstencil+lambdamax(i)*vsv)
            vr(0)=vstencil(0)*cr(0,0)+vstencil(1)*cr(1,0)
     $           +vstencil(2)*cr(2,0)
            vr(1)=vstencil(-1)*cr(0,1)+vstencil(0)*cr(1,1)
     $           +vstencil(1)*cr(2,1)
            vr(2)=vstencil(-2)*cr(0,2)+vstencil(-1)*cr(1,2)
     $           +vstencil(0)*cr(2,2)
            belta(0)=13.d0/12.*(vstencil(0)-2.*vstencil(1)
     $              +vstencil(2))**2+0.25*(3.*vstencil(0)
     $              -4.*vstencil(1)+vstencil(2))**2
            belta(1)=13.d0/12.*(vstencil(-1)-2.*vstencil(0)
     $              +vstencil(1))**2+0.25*(vstencil(-1)
     $              -vstencil(1))**2
            belta(2)=13.d0/12.*(vstencil(-2)-2.*vstencil(-1)
     $              +vstencil(0))**2+0.25*(vstencil(-2)
     $              -4.*vstencil(-1)+3.*vstencil(0))**2
            alpha(0)=dr(1,0)/(1.0e-7+belta(0))**2
            alpha(1)=dr(1,1)/(1.0e-7+belta(1))**2
            alpha(2)=dr(1,2)/(1.0e-7+belta(2))**2
            omiga(0)=alpha(0)/(alpha(0)+alpha(1)+alpha(2))
            omiga(1)=alpha(1)/(alpha(0)+alpha(1)+alpha(2))
            omiga(2)=1.0d0-omiga(0)-omiga(1)
            vflux(i)=omiga(0)*vr(0)+omiga(1)*vr(1)+omiga(2)*vr(2)
c     -------- caculate the g(i+1/2)- ------------
            vstencil=g(i,-1:3)
            vsv=v(i,-1:3)

            vstencil=0.5*(vstencil-lambdamax(i)*vsv)
            vr(0)=vstencil(0)*cr(0,-1)+vstencil(1)*cr(1,-1)
     $           +vstencil(2)*cr(2,-1)
            vr(1)=vstencil(-1)*cr(0,0)+vstencil(0)*cr(1,0)
     $           +vstencil(1)*cr(2,0)
            vr(2)=vstencil(-2)*cr(0,1)+vstencil(-1)*cr(1,1)
     $           +vstencil(0)*cr(2,1)
            belta(0)=13.d0/12.*(vstencil(0)-2.*vstencil(1)
     $              +vstencil(2))**2+0.25*(3.*vstencil(0)
     $              -4.*vstencil(1)+vstencil(2))**2
            belta(1)=13.d0/12.*(vstencil(-1)-2.*vstencil(0)
     $              +vstencil(1))**2+0.25*(vstencil(-1)
     $              -vstencil(1))**2
            belta(2)=13.d0/12.*(vstencil(-2)-2.*vstencil(-1)
     $              +vstencil(0))**2+0.25*(vstencil(-2)
     $              -4.*vstencil(-1)+3.*vstencil(0))**2
            alpha(0)=dr(1,0)/(1.0e-7+belta(0))**2
            alpha(1)=dr(1,1)/(1.0e-7+belta(1))**2
            alpha(2)=dr(1,2)/(1.0e-7+belta(2))**2
            omiga(0)=alpha(0)/(alpha(0)+alpha(1)+alpha(2))
            omiga(1)=alpha(1)/(alpha(0)+alpha(1)+alpha(2))
            omiga(2)=1.0d0-omiga(0)-omiga(1)
            alpha(0)=omiga(0)*vr(0)+omiga(1)*vr(1)+omiga(2)*vr(2)
            vflux(i)=vflux(i)+alpha(0)
        enddo        
        flux(:,ix)=matmul(rightv,vflux)        
      enddo 
      end
c     =====================================================

      subroutine output(u,x,N)
      integer N,i 
      real(8) pi,gama,cv
      real(8) u(3,-3:N+3),x(-3:N+3)
      real(8) ro(0:N),v(0:N),p(0:N),T(0:N),u0(0:N)
      common /group2/ gama,cv
      pi=acos(-1.d0)

       open(unit=02,file='./WENO_Sod.dat',status='unknown')
       write(02,101)
       write(02,102)
       write(02,103) N+1

      do i=0,N
        ro(i)=u(1,i)
        p(i)=(gama-1.d0)*(u(3,i)-0.5*u(2,i)*u(2,i)/u(1,i))
        v(i) =u(2,i)/u(1,i)
        ! T(i)=p(i)/((gama-1.d0)*ro(i))
        write(02,100) x(i),ro(i),p(i),v(i) !,T(i)
      enddo

100    format(2x,10(e12.6,'      '))
101    format('Title="Sod Shock Tube(WENO Scheme)"')
102    format('Variables=x,rho,p,u')
103    format('zone',1x,'i=',1x,i5,2x,'f=point')

       close(02)
       write(*,*) 'Data export to ./WENO_Sod.dat file!'

      end

c     =====================================================
      subroutine ghostp(u,N,dx)
      integer N
      real(8) u(3,-3:N+3),dx
      real(8) leftv(3,3),vch(3,5),vbd(3,0:4),ubd(3,0:4),lambda(3)
      real(8) ro,vel,p,e,gama,cv,c
      real(8) beltabd(0:4),alphabd(0:4),omigabd(0:4)
      real(8) u0,u1,u2,u3,u4,u5,dbd(0:4),mu,zi
      real(8) dkprdxk(0:4,0:4),rigv(3,3)
      common /group2/ gama,cv

c     caculate the ghost points
c     ghost points u(n+1),u(n+2),u(n+3)
      ro=u(1,N)
      vel=u(2,N)/ro
      p=(gama-1.d0)*(u(3,N)-0.5*ro*vel*vel)
      e=u(3,N)/ro
      c=sqrt(gama*p/ro)
      lambda(1)=vel-c
      lambda(2)=vel
      lambda(3)=vel+c
c     right eigenvector
      rigv=1.d0
      rigv(2,1)=lambda(1);rigv(3,1)=vel*vel*0.5-vel*c+c*c/(gama-1.d0)
      rigv(2,2)=vel;rigv(3,2)=vel*vel*0.5
      rigv(2,3)=lambda(3);rigv(3,3)=vel*vel*0.5+vel*c+c*c/(gama-1.d0)
      call leftvector(leftv,rigv)

      inflownum=0
      if(lambda(1)<0) inflownum=inflownum+1
      if(lambda(2)<0) inflownum=inflownum+1
      if(lambda(3)<0) inflownum=inflownum+1
      vch=matmul(leftv,u(:,N-4:N))

      do i=1,3
        u0=vch(i,5);u1=vch(i,4);u2=vch(i,3);u3=vch(i,2);u4=vch(i,1)
        beltabd(0)=dx*dx
        beltabd(1)=(u0-u1)*(u0-u1)
        beltabd(2)=(61.d0*u0*u0-196.d0*u0*u1+74.d0*u0*u2+
     $              160.d0*u1*u1-124.d0*u1*u2+25*u2*u2)/12.d0
        beltabd(3)=39.63333333333333d0*u1*u3-144.9833333333333d0*u1*u2
     $            -32.38333333333333d0*u2*u3-u0*(77.71666666666667d0*u1
     $            -60.3d0*u2+16.29444444444444d0*u3)+16.85555555555556d0
     $            *u0*u0+91.53333333333333d0*u1*u1+58.53333333333333d0
     $            *u2*u2+4.522222222222222d0*u3*u3
        beltabd(4)=47.17781084656085d0*u0*u0-290.0627645502646d0*u0*u1
     $            +345.4406746031746d0*u0*u2-190.1155423280423d0*u0*u3 
     $            +40.38201058201058d0*u0*u4+455.7616402116402d0*u1*u1
     $            -1104.154365079365d0*u1*u2+613.9899470899471d0*u1*u3 
     $            -131.2960978835979d0*u1*u4+677.7761904761905d0*u2*u2 
     $            -760.4043650793651d0*u2*u3+163.5656746031746d0*u2*u4
     $            +214.6505291005291d0*u3*u3-92.77109788359788d0*u3*u4 
     $            +10.05975529100529d0*u4*u4
        dbd(0)=dx**4;dbd(1)=dx**3;dbd(2)=dx*dx;dbd(3)=dx
        dbd(4)=1-(dbd(0)+dbd(1)+dbd(2)+dbd(3))
        do j=0,4
            alphabd(j)=dbd(j)/(1.e-6+beltabd(j))**3
        enddo
        do j=0,4
            omigabd(j)=alphabd(j)/(alphabd(0)+alphabd(1)
     $                +alphabd(2)+alphabd(3)+alphabd(4))
        enddo
        dkprdxk=0.d0
        dkprdxk(0,0)=u0
        dkprdxk(1,0)=1.5*u0 - 0.5*u1
        dkprdxk(2,0)=1.875*u0 - 1.25*u1 + 0.375*u2
        dkprdxk(3,0)=2.1875*u0 - 2.1875*u1 + 1.3125*u2 - 0.3125*u3
        dkprdxk(4,0)=2.4609375*u0-3.28125*u1+2.953125*u2-1.40625*u3
     $              +0.2734375*u4
        dkprdxk(1,1)=u0/dx-u1/dx
        dkprdxk(2,1)=2.0*u0/dx - 3.0*u1/dx + u2/dx
        dkprdxk(3,1)=2.958333333333333d0*u0/dx - 5.875*u1/dx 
     $              +3.875*u2/dx - 0.9583333333333333d0*u3/dx
        dkprdxk(4,1)=3.875*u0/dx - 9.541666666666667d0*u1/dx 
     $              +9.375*u2/dx - 4.625*u3/dx 
     $              +0.9166666666666667d0*u4/dx
        dkprdxk(2,2)=u0/(dx*dx) - 2.0*u1/(dx*dx) + u2/(dx*dx)
        dkprdxk(3,2)=2.5*u0/(dx*dx) - 6.5*u1/(dx*dx) 
     $              +5.5*u2/(dx*dx) - 1.5*u3/(dx*dx)
        dkprdxk(4,2)=4.291666666666667d0*u0/(dx*dx) 
     $              -13.66666666666667d0*u1/(dx*dx) + 16.25*u2/(dx*dx) 
     $              -8.666666666666667d0*u3/(dx*dx) 
     $              +1.791666666666667d0*u4/(dx*dx)
        dkprdxk(3,3)=u0/(dx*dx*dx) - 3.0*u1/(dx*dx*dx) 
     $              +3.0*u2/(dx*dx*dx) - u3/(dx*dx*dx)
        dkprdxk(4,3)=3.0*u0/(dx*dx*dx) - 11.0*u1/(dx*dx*dx) 
     $              +15.0*u2/(dx*dx*dx)-9.0*u3/(dx*dx*dx)
     $              +2.0*u4/(dx*dx*dx)
        dkprdxk(4,4)=u0/(dx*dx*dx*dx) - 4.0*u1/(dx*dx*dx*dx) 
     $              +6.0*u2/(dx*dx*dx*dx) - 4.0*u3/(dx*dx*dx*dx)
     $              + u4/(dx*dx*dx*dx)
        vbd(i,:)=matmul(omigabd,dkprdxk)
      enddo

      if(inflownum/=1)then
c        write(*,*) 'WARNING!'
c        write(*,*) lambda
      endif

      ubd(2,0)=0.d0
      ubd(1,0)=(vbd(2,0)*leftv(3,3)-vbd(3,0)*leftv(2,3))
     $        /(leftv(2,1)*leftv(3,3)-leftv(2,3)*leftv(3,1))
      ubd(3,0)=(vbd(3,0)*leftv(2,1)-vbd(2,0)*leftv(3,1))
     $        /(leftv(2,1)*leftv(3,3)-leftv(2,3)*leftv(3,1))
      mu=0.5*(gama-3.)*(ubd(2,0)/ubd(1,0))**2*leftv(2,2)*leftv(3,3)
     $  +(3.-gama)*ubd(2,0)/ubd(1,0)*leftv(2,3)*leftv(3,1)
     $  +leftv(2,1)*leftv(3,2)*(gama-1.)-(gama-1.)*leftv(2,2)*leftv(3,1)
     $  -0.5*(gama-3.)*(ubd(2,0)/ubd(1,0))**2*leftv(2,3)*leftv(3,2)
     $  -(3.-gama)*ubd(2,0)/ubd(1,0)*leftv(2,1)*leftv(3,3)
      zi=(3.-gama)*ubd(2,0)/ubd(1,0)*leftv(2,3)*vbd(3,1)
     $  +vbd(2,1)*leftv(3,2)*(gama-1.)
     $  -(gama-1.)*leftv(2,2)*vbd(3,1)
     $  -(3.-gama)*ubd(2,0)/ubd(1,0)*vbd(2,1)*leftv(3,3)
      ubd(1,1)=zi/mu
      zi=0.5*(gama-3.)*(ubd(2,0)/ubd(1,0))**2*vbd(2,1)*leftv(3,3)
     $  +leftv(2,1)*vbd(3,1)*(gama-1.)-(gama-1.)*vbd(2,1)*leftv(3,1)
     $  -0.5*(gama-3.)*(ubd(2,0)/ubd(1,0))**2*leftv(2,3)*vbd(3,1)
      ubd(2,1)=zi/mu
      zi=0.5*(gama-3.)*(ubd(2,0)/ubd(1,0))**2*leftv(2,2)*vbd(3,1)
     $  +(3.-gama)*ubd(2,0)/ubd(1,0)*vbd(2,1)*leftv(3,1)
     $  -0.5*(gama-3.)*(ubd(2,0)/ubd(1,0))**2*vbd(2,1)*leftv(3,2)
     $  -(3.-gama)*ubd(2,0)/ubd(1,0)*leftv(2,1)*vbd(3,1)
      ubd(3,1)=zi/mu
      
      ubd(:,2:4)=matmul(rigv,vbd(:,2:4))
      do i=1,3
      mu=(0.5-i)*dx
      u(:,N+i)=ubd(:,0)+mu*ubd(:,1)+0.5*mu*mu*ubd(:,2)
     $           +mu*mu*mu/6.*ubd(:,3)+mu**4/24.*ubd(:,4)
      enddo
c     ----------------------------------------------------------
c     ghost points u(-1),u(-2),u(-3)
      ro=u(1,0)
      vel=u(2,0)/ro
      p=(gama-1.d0)*(u(3,0)-0.5*ro*vel*vel)
      e=u(3,0)/ro
      c=sqrt(gama*p/ro)
      lambda(1)=vel-c
      lambda(2)=vel
      lambda(3)=vel+c
c     right eigenvector
      rigv=1.d0
      rigv(2,1)=lambda(1)
      rigv(3,1)=vel*vel*0.5-vel*c+c*c/(gama-1.d0)
      rigv(2,2)=vel
      rigv(3,2)=vel*vel*0.5
      rigv(2,3)=lambda(3)
      rigv(3,3)=vel*vel*0.5+vel*c+c*c/(gama-1.d0)
      call leftvector(leftv,rigv)

      inflownum=0
      if(lambda(1)>0) inflownum=inflownum+1
      if(lambda(2)>0) inflownum=inflownum+1
      if(lambda(3)>0) inflownum=inflownum+1
      vch=matmul(leftv,u(:,0:4))
      
      do i=1,3
        u0=vch(i,1);u1=vch(i,2);u2=vch(i,3);u3=vch(i,4);u4=vch(i,5)
        beltabd(0)=dx*dx
        beltabd(1)=(u0-u1)*(u0-u1)
        beltabd(2)=(61.d0*u0*u0-196.d0*u0*u1+74.d0*u0*u2+
     $              160.d0*u1*u1-124.d0*u1*u2+25*u2*u2)/12.d0
        beltabd(3)=39.63333333333333d0*u1*u3-144.9833333333333d0*u1*u2
     $            -32.38333333333333d0*u2*u3-u0*(77.71666666666667d0*u1
     $            -60.3d0*u2+16.29444444444444d0*u3)+16.85555555555556d0
     $            *u0*u0+91.53333333333333d0*u1*u1+58.53333333333333d0
     $            *u2*u2+4.522222222222222d0*u3*u3
        beltabd(4)=47.17781084656085d0*u0*u0-290.0627645502646d0*u0*u1
     $            +345.4406746031746d0*u0*u2-190.1155423280423d0*u0*u3 
     $            +40.38201058201058d0*u0*u4+455.7616402116402d0*u1*u1
     $            -1104.154365079365d0*u1*u2+613.9899470899471d0*u1*u3 
     $            -131.2960978835979d0*u1*u4+677.7761904761905d0*u2*u2 
     $            -760.4043650793651d0*u2*u3+163.5656746031746d0*u2*u4
     $            +214.6505291005291d0*u3*u3-92.77109788359788d0*u3*u4 
     $            +10.05975529100529d0*u4*u4
        dbd(0)=dx**4;dbd(1)=dx**3;dbd(2)=dx*dx;dbd(3)=dx
        dbd(4)=1-(dbd(0)+dbd(1)+dbd(2)+dbd(3))
        do j=0,4
            alphabd(j)=dbd(j)/(1.e-6+beltabd(j))**3
        enddo
        do j=0,4
            omigabd(j)=alphabd(j)/(alphabd(0)+alphabd(1)
     $                +alphabd(2)+alphabd(3)+alphabd(4))
        enddo
        dkprdxk=0.d0
        dkprdxk(0,0)=u0
        dkprdxk(1,0)=1.5*u0 - 0.5*u1
        dkprdxk(2,0)=1.875*u0 - 1.25*u1 + 0.375*u2
        dkprdxk(3,0)=2.1875*u0 - 2.1875*u1 + 1.3125*u2 - 0.3125*u3
        dkprdxk(4,0)=2.4609375*u0-3.28125*u1+2.953125*u2-1.40625*u3
     $              +0.2734375*u4
        dkprdxk(1,1)=u1/dx-u0/dx
        dkprdxk(2,1)=3.0*u1/dx - 2.0*u0/dx - u2/dx
        dkprdxk(3,1)=-2.958333333333333d0*u0/dx + 5.875*u1/dx 
     $              -3.875*u2/dx + 0.9583333333333333d0*u3/dx
        dkprdxk(4,1)=-3.875*u0/dx + 9.541666666666667d0*u1/dx 
     $              -9.375*u2/dx + 4.625*u3/dx 
     $              -0.9166666666666667d0*u4/dx
        dkprdxk(2,2)=u0/(dx*dx) - 2.0*u1/(dx*dx) + u2/(dx*dx)
        dkprdxk(3,2)=2.5*u0/(dx*dx) - 6.5*u1/(dx*dx) 
     $              +5.5*u2/(dx*dx) - 1.5*u3/(dx*dx)
        dkprdxk(4,2)=4.291666666666667d0*u0/(dx*dx) 
     $              -13.66666666666667d0*u1/(dx*dx) + 16.25*u2/(dx*dx) 
     $              -8.666666666666667d0*u3/(dx*dx) 
     $              +1.791666666666667d0*u4/(dx*dx)
        dkprdxk(3,3)=-u0/(dx*dx*dx) + 3.0*u1/(dx*dx*dx) 
     $              -3.0*u2/(dx*dx*dx) + u3/(dx*dx*dx)
        dkprdxk(4,3)=-3.0*u0/(dx*dx*dx) + 11.0*u1/(dx*dx*dx) 
     $              -15.0*u2/(dx*dx*dx)+9.0*u3/(dx*dx*dx)
     $              -2.0*u4/(dx*dx*dx)
        dkprdxk(4,4)=u0/(dx*dx*dx*dx) - 4.0*u1/(dx*dx*dx*dx) 
     $              +6.0*u2/(dx*dx*dx*dx) - 4.0*u3/(dx*dx*dx*dx)
     $              + u4/(dx*dx*dx*dx)
        vbd(i,:)=matmul(omigabd,dkprdxk)
      enddo

      if(inflownum/=1)then
c        write(*,*) 'WARNING!'
c        write(*,*) lambda
      endif
c     lambda3>0 as the inflow variable
      ubd(2,0)=0.d0
      ubd(1,0)=(vbd(1,0)*leftv(2,3)-vbd(2,0)*leftv(1,3))
     $        /(leftv(1,1)*leftv(2,3)-leftv(1,3)*leftv(2,1))
      ubd(3,0)=(vbd(2,0)*leftv(1,1)-vbd(1,0)*leftv(2,1))
     $        /(leftv(1,1)*leftv(2,3)-leftv(1,3)*leftv(2,1))

      mu=-0.5*(gama-3.)*(ubd(2,0)/ubd(1,0))**2*leftv(2,2)*leftv(1,3)
     $  -(3.-gama)*ubd(2,0)/ubd(1,0)*leftv(2,3)*leftv(1,1)
     $  -leftv(2,1)*leftv(1,2)*(gama-1.)
     $  +(gama-1.)*leftv(2,2)*leftv(1,1)
     $  +0.5*(gama-3.)*(ubd(2,0)/ubd(1,0))**2*leftv(2,3)*leftv(1,2)
     $  +(3.-gama)*ubd(2,0)/ubd(1,0)*leftv(2,1)*leftv(1,3)

      zi=(3.-gama)*ubd(2,0)/ubd(1,0)*leftv(1,3)*vbd(2,1)
     $  +vbd(1,1)*leftv(2,2)*(gama-1.)
     $  -(gama-1.)*leftv(1,2)*vbd(2,1)
     $  -(3.-gama)*ubd(2,0)/ubd(1,0)*vbd(1,1)*leftv(2,3)
      ubd(1,1)=zi/mu

      zi=0.5*(gama-3.)*(ubd(2,0)/ubd(1,0))**2*vbd(1,1)*leftv(2,3)
     $  +leftv(1,1)*vbd(2,1)*(gama-1.)
     $  -(gama-1.)*vbd(1,1)*leftv(2,1)
     $  -0.5*(gama-3.)*(ubd(2,0)/ubd(1,0))**2*leftv(1,3)*vbd(2,1)
      ubd(2,1)=zi/mu

      zi=0.5*(gama-3.)*(ubd(2,0)/ubd(1,0))**2*leftv(1,2)*vbd(2,1)
     $  +(3.-gama)*ubd(2,0)/ubd(1,0)*vbd(1,1)*leftv(2,1)
     $  -0.5*(gama-3.)*(ubd(2,0)/ubd(1,0))**2*vbd(1,1)*leftv(2,2)
     $  -(3.-gama)*ubd(2,0)/ubd(1,0)*leftv(1,1)*vbd(2,1)
      ubd(3,1)=zi/mu
      
      ubd(:,2:4)=matmul(rigv,vbd(:,2:4))
      do i=-3,-1
      mu=(0.5+i)*dx
      u(:,i)=ubd(:,0)+mu*ubd(:,1)+0.5*mu*mu*ubd(:,2)
     $           +mu*mu*mu/6.*ubd(:,3)+mu**4/24.*ubd(:,4)
      enddo

c     use the code below as periodic BCs
c      u(:,0)=u(:,N)
c      u(:,-3:-1)=u(:,N-3:N-1)
c      u(:,N+1:N+3)=u(:,1:3)

      end
c     =====================================================

c     =====================================================
      subroutine leftvector(rightv,leftv)
      real(8) leftv(3,3),rightv(3,3)
      real(8) mu

      mu=leftv(1,1)*leftv(2,2)*leftv(3,3)
     $  +leftv(1,2)*leftv(2,3)*leftv(3,1)
     $  +leftv(2,1)*leftv(3,2)*leftv(1,3)
     $  -leftv(1,3)*leftv(2,2)*leftv(3,1)
     $  -leftv(1,1)*leftv(2,3)*leftv(3,2)
     $  -leftv(1,2)*leftv(2,1)*leftv(3,3)
      rightv(1,1)= (leftv(2,2)*leftv(3,3)-leftv(2,3)*leftv(3,2))/mu
      rightv(2,2)= (leftv(1,1)*leftv(3,3)-leftv(1,3)*leftv(3,1))/mu
      rightv(3,3)= (leftv(2,2)*leftv(1,1)-leftv(2,1)*leftv(1,2))/mu
      rightv(2,1)=-(leftv(2,1)*leftv(3,3)-leftv(2,3)*leftv(3,1))/mu
      rightv(3,1)= (leftv(2,1)*leftv(3,2)-leftv(2,2)*leftv(3,1))/mu
      rightv(1,2)=-(leftv(1,2)*leftv(3,3)-leftv(1,3)*leftv(3,2))/mu
      rightv(3,2)=-(leftv(1,1)*leftv(3,2)-leftv(1,2)*leftv(3,1))/mu
      rightv(1,3)= (leftv(1,2)*leftv(2,3)-leftv(2,2)*leftv(1,3))/mu
      rightv(2,3)=-(leftv(1,1)*leftv(2,3)-leftv(2,1)*leftv(1,3))/mu  
      end
c     =====================================================

c     =====================================================
