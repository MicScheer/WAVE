*CMZ :  3.05/23 12/11/2018  16.49.45  by  Michael Scheer
*-- Author :    Michael Scheer   12/11/2018
      subroutine mrad_fringe_quintic_spline_msh(
     &  xin,y,z,bx,by,bz,ax,ay,az,fint,gap,fringe,istatus)

      implicit none

      double precision, intent(in) :: xin,y,z,gap,fint
      double precision, intent(out) :: bx,by,bz,fringe,ax,ay,az !ax,ay,az dummy
      double precision fb,fa,fc,fringe2,fringe3,fringe4,fringe5,
     &  x2,x3,x4,x5,x,y2,y3

      integer, intent(out) :: istatus

      include 'phyconparam.cmn'

      istatus=0

      x=xin

      bx=0.0d0
      by=0.0d0
      bz=0.0d0

      ax=0.0d0
      ay=0.0d0
      az=0.0d0

      if (gap.le.0.0d0.or.fint.le.0.0d0) then
        istatus=-1
        goto 9999
      endif

      fringe=231.0d0*fint*gap/25.0d0

      if(x.lt.0.0d0) then
        return
      else if (x.gt.fringe) then
        x=fringe
      endif

      fringe2=fringe*fringe
      fringe3=fringe2*fringe
      fringe4=fringe2*fringe2
      fringe5=fringe3*fringe2

      x2=x*x
      x3=x2*x
      x4=x2*x2
      x5=x3*x2
      y2=y*y
      y3=y2*y

      fa=10.0d0/fringe3
      fb=-15.0d0/fringe4
      fc=6.0d0/fringe5

      bx=y*(3.0d0*Fa*x2+4.0d0*fb*x3+5.0d0*fc*x4)
c     &  +y3*(-fa-4.0d0*fb*x-10.0d0*fc*x2) ! it is not Maxwell conform
c      bx=x2*y*(-3.0d0*fa-4.0d0*fb*x-5.0d0*fc*x2)
      by=(fa*x3+fb*x4+fc*x5)-y2*x*(3.0d0*fa+6.0d0*fb*x+10.0d0*fc*x2)

      bx=y*(3.0d0*Fa*x2+4.0d0*fb*x3+5.0d0*fc*x4)
     &  +y3*(-fa-4.0d0*fb*x-10.0d0*fc*x2)
      by=(fa*x3+fb*x4+fc*x5)-y2*x*(3.0d0*fa+6.0d0*fb*x+10.0d0*fc*x2)
     &  +y**4*(fb+5.0d0*fc*x) !this makes it Maxwell conform

9999  if (istatus.ne.0) then
        print*,""
        print*,"*** Error in mrad_quintic_fringe_spline_msh:",istatus
        print*,""
      endif

      return
      end
