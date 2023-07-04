*CMZ :  4.00/11 21/04/2021  12.07.45  by  Michael Scheer
*-- Author :    Michael Scheer   19/04/2021
      subroutine trajectory_to_bfield(nstep,rmass,charge,ebeam,xyz,bxyz,istat)

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      double precision, dimension (:),  allocatable ::
     &  x,y,z,xp,yp,zp,xpp,ypp,zpp,aa,bb,cc,c,t

      double precision rmass,dmass,charge,ebeam,xyz(3,nstep),bxyz(3,nstep),
     &  beta,gamma,a,b,ax,ay,az,bn,v0,bx,by,bz

      integer nstep,istat,i

      ! x(i) = xyz(1,i)
      ! y(i) = xyz(2,i)
      ! z(i) = xyz(3,i)

      istat=0

      gamma=ebeam/emassg1
      beta=dsqrt((1.0d0-1.0d0/gamma)*(1.0d0+1.0d0/gamma))
      v0=clight1*beta
      dmass=rmass*gamma

      allocate(t(nstep),
     &  x(nstep),y(nstep),z(nstep),
     &  xp(nstep),yp(nstep),zp(nstep),
     &  xpp(nstep),ypp(nstep),zpp(nstep),
     &  aa(nstep),bb(nstep),cc(nstep),c(nstep))

      x=xyz(1,:)
      y=xyz(2,:)
      z=xyz(3,:)

      t(1)=0.0d0
      do i=2,nstep
        t(i)=t(i-1)+sqrt((x(i)-x(i-1))**2+(y(i)-y(i-1))**2+(z(i)-z(i-1))**2)/v0
      enddo

      call util_spline_coef_deriv(t,x,nstep,9999.0d0,9999.0d0,xp,xpp,
     &  aa,bb,cc,c)
      call util_spline_coef_deriv(t,y,nstep,9999.0d0,9999.0d0,yp,ypp,
     &  aa,bb,cc,c)
      call util_spline_coef_deriv(t,z,nstep,9999.0d0,9999.0d0,zp,zpp,
     &  aa,bb,cc,c)

      do i=1,nstep
        ax=xpp(i)
        ay=ypp(i)
        az=zpp(i)
        a=sqrt(ax*ax+ay*ay+az*az)
        b=a/v0*dmass/charge
        bx=(yp(i)*zpp(i)-zp(i)*ypp(i))
        by=(zp(i)*xpp(i)-xp(i)*zpp(i))
        bz=(xp(i)*ypp(i)-yp(i)*xpp(i))
        bn=b/sqrt(bx*bx+by*by+bz*bz)
        bxyz(1,i)=bx*bn
        bxyz(2,i)=by*bn
        bxyz(3,i)=bz*bn
      enddo

      deallocate(t,x,y,z,xp,yp,zp,xpp,ypp,zpp,aa,bb,cc,c)

      return
      end
