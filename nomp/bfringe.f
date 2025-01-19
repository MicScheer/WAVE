*CMZ :  4.01/07 13/01/2025  22.09.50  by  Michael Scheer
*-- Author :    Michael Scheer   23/11/2024
      subroutine bfringe(cbmodel,x,y,z,b0,bx,by,bz,fringe,fa,fb,fc,icharge,istatus)

      implicit none

cdbuff dcbuff

      double precision, intent(in) :: x,y,z,b0,fringe,fa,fb,fc
      double precision, intent(out) :: bx,by,bz
      double precision x2,x3,x4,x5,y2,y3

      character(32), intent(in) :: cbmodel

      integer, intent(in) :: icharge
      integer, intent(out) :: istatus

      include 'phyconparam.cmn'

      istatus=0

      bx=0.0d0
      by=0.0d0
      bz=0.0d0

      if(x.lt.0.0d0) return

      if (cbmodel.eq."quintic-spline") then
        x2=x*x
        x3=x2*x
        x4=x2*x2
        x5=x3*x2
        y2=y*y
        y3=y2*y
        bx=y*(3.0d0*Fa*x2+4.0d0*fb*x3+5.0d0*fc*x4)
     &    +y3*(-fa-4.0d0*fb*x-10.0d0*fc*x2) !This term is not Maxwell conform
c             by2=(fa*x3+fb*x4+fc*x5)+y2*x*(3.0d0*fa+6.0d0*fb*x+10.0d0*fc*x2) ! The sign seems to be wrong in the manual
        by=(fa*x3+fb*x4+fc*x5)-y2*x*(3.0d0*fa+6.0d0*fb*x+10.0d0*fc*x2)
      else if (cbmodel.eq.'cubic-spline') then
        x2=x*x
        x3=x2*x
        bx=(2.0d0*fa*x+3.0d0*fb*x2)*y
        by=(fa*x2+fb*x3+y**2*(-fa-3.0d0*fb*x))
      else if (cbmodel.eq.'linear') then
        bx=y*fa
        by=x*fa
      else
        stop "*** Bad model in bfringe ***"
      endif

      bx=bx*b0
      by=by*b0

      return
      end
