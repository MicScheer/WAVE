*CMZ :  3.05/26 06/12/2018  11.57.49  by  Michael Scheer
*CMZ :  3.05/23 27/11/2018  11.49.31  by  Michael Scheer
*-- Author :    Michael Scheer   12/11/2018
      subroutine mrad_fringe_linear(
     &  xin,y,z,bx,by,bz,ax,ay,az,fint,gap,fringe,istatus)

      implicit none

      double precision, intent(in) :: xin,y,z,gap,fint
      double precision, intent(out) :: bx,by,bz,fringe,ax,ay,az !ax,ay,az dummy
      double precision fa,x

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

      fringe=6.0d0*fint*gap

      if(x.lt.0.0d0) then
        by=0.0d0
        return
      else if (x.gt.fringe) then
        by=1.0d0
      else
        fa=1.0d0/fringe
        bx=y*fa
        by=x*fa
      endif

9999  continue

      return
      end
