*CMZ :  3.06/00 18/01/2019  12.48.01  by  Michael Scheer
*CMZ :  3.05/28 04/01/2019  17.40.38  by  Michael Scheer
*CMZ :  3.05/26 06/12/2018  11.57.05  by  Michael Scheer
*CMZ :  3.05/23 12/11/2018  13.44.45  by  Michael Scheer
*-- Author :    Michael Scheer   12/11/2018
      subroutine mrad_fringe_cubic_spline(
     &  xin,y,z,bx,by,bz,ax,ay,az,fint,gap,fringe,istatus)

      implicit none

      double precision, intent(in) :: xin,y,z,gap,fint
      double precision, intent(out) :: bx,by,bz,fringe,ax,ay,az !ax,ay,az dummy
      double precision fb,fa,fringe2,fringe3,x2,x3,x

      integer, intent(out) :: istatus

      include 'phyconparam.cmn'

      istatus=0

      x=xin

      bx=0.0d0
      by=1.0d0
      bz=0.0d0

      ax=0.0d0
      ay=0.0d0
      az=0.0d0

      if (gap.le.0.0d0.or.fint.le.0.0d0) then
        goto 9999
      endif

      !Int(y=0,by,0->fringe)=1/2

      fringe=70.0d0/9.0d0*fint*gap

      if(x.lt.0.0d0) then
        by=0.0d0
        return
      else if (x.gt.fringe) then
        by=1.0d0
      else

        fringe2=fringe*fringe
        fringe3=fringe2*fringe

        x2=x*x
        x3=x2*x

        fb=-2.0d0/fringe3
        fa=3.0d0/fringe2

        bx=(2.0d0*fa*x+3.0d0*fb*x2)*y
        by=fa*x2+fb*x3+y**2*(-fa-3.0d0*fb*x)

      endif

9999  continue

      return
      end
