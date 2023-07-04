*CMZ :  3.05/23 21/11/2018  13.03.11  by  Michael Scheer
*-- Author :    Michael Scheer   21/11/2018
      subroutine util_atan(x,y,grad,rad)

      implicit none

      double precision x,y,grad,rad

      rad=atan2(y,x)
      if (rad.lt.0.0d0) rad=rad+6.283185307180d0
      grad=rad*57.295779513082323d0

      return
      end
