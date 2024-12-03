*CMZ : 00.00/16 08/07/2014  15.30.50  by  Michael Scheer
*CMZ : 00.00/15 18/01/2013  16.42.57  by  Michael Scheer
*-- Author :    Michael Scheer   18/01/2013
      subroutine util_circle(x,y,x0,y0,r)

      implicit none

      double precision x(3),y(3),x0,y0,r

      r=(2.0d0*(x(1)*y(2)-x(1)*y(3)-x(2)*y(1)+x(2)*y(3)+x(3)*y(1)-x(3)*y(2)))
      x0=0.0d0
      y0=0.0d0

      if (abs(r).gt.1.0d-10) then
        x0=(x(1)**2*y(2)-x(1)**2*y(3)-x(2)**2*y(1)+x(2)**2*y(3)+
     &    x(3)**2*y(1)-x(3)**2*y(2)+y(1)**2*y(2)-y(1)**2*y(3)-y(1)*y(2)**2+y(1)*y(3)**2+
     &    y(2)**2*y(3)-y(2)*y(3)**2)/r
      else
        r=0.0d0
        return
      endif

      y0=(-x(1)**2*x(2)+x(1)**2*x(3)+x(1)*x(2)**2-x(1)*x(3)**2+x(1)*y(2)**2-x(1)*
     &  y(3)**2-x(2)**2*x(3)+x(2)*x(3)**2-x(2)*y(1)**2+x(2)*y(3)**2+
     &  x(3)*y(1)**2-x(3)*y(2)**2)/r

      r=sqrt((x(1)-x0)**2+(y(1)-y0)**2)

      return
      end
