*CMZ :          16/06/2017  16.47.35  by  Michael Scheer
*-- Author :    Michael Scheer   16/06/2017

      subroutine util_weed_line(n,x,y,z,n1,tiny)

c+PATCH,//UTIL/FOR
c+DECK,util_weed_line.

      implicit none

      integer n,n1

      double precision x(n),y(n),z(n),xmin,ymin,zmin,xmax,ymax,zmax,
     &  v(3),vn,tiny

      double precision, dimension (:), allocatable :: xb,yb,zb

      integer i

      stop "NOT READY "

      n1=n

      if (n.lt.3) then
        return
      endif

      allocate(xb(n),yb(n),zb(n))

      xmax=-1.0d30
      xmin=1.0d30
      ymax=-1.0d30
      ymin=1.0d30
      zmax=-1.0d30
      zmin=1.0d30

      do i=1,n
        if (x(i).lt.xmin) xmin=x(i)
        if (x(i).gt.xmax) xmax=x(i)
        if (y(i).lt.ymin) ymin=y(i)
        if (y(i).gt.ymax) ymax=y(i)
        if (z(i).lt.zmin) zmin=z(i)
        if (z(i).gt.zmax) zmax=z(i)
      enddo

      if (xmax.ne.xmin) then
        xb=(x-xmin)/(xmax-xmin)
      else
        xb=0.0d0
      endif

      if (ymax.ne.ymin) then
        yb=(y-ymin)/(ymax-ymin)
      else
        yb=0.0d0
      endif

      do i=2,n
        v(1)=x(i)-x(n)
        v(2)=y(i)-x(n)
        v(3)=z(i)-x(n)
      enddo

9999  deallocate(xb,yb,zb)

      return
      end
