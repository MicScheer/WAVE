*CMZ :          04/06/2019  11.56.47  by  Michael Scheer
*-- Author :    Michael Scheer   12/04/2016
      subroutine util_sort_triangle(n,x,y,itri)

      implicit none

      double precision x(n),y(n),
     &  x1,y1,x2,y2,x3,y3,dx1,dy1,dx2,dy2,dz

      integer n,itri(3),i

      x1=x(itri(1))
      y1=y(itri(1))
      x2=x(itri(2))
      y2=y(itri(2))
      x3=x(itri(3))
      y3=y(itri(3))

      dx1=x2-x1
      dy1=y2-y1
      dx2=x3-x2
      dy2=y3-y2

      dz=dx1*dy2-dx2*dy1

      if (dz.lt.0.0d0) then
        i=itri(3)
        itri(3)=itri(2)
        itri(2)=i
      endif

      return
      end
