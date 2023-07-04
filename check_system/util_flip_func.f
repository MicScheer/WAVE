*CMZ : 00.00/11 29/03/2011  15.49.35  by  Michael Scheer
*CMZ : 00.00/00 11/01/95  11.41.04  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_flip_func(n,x,y,wx,wy)

c flips (x1,y1),...,(xn,yn) to (xn,yn),...,(x1,y1)

      implicit none

      integer n,i,in
      double precision x(n),y(n),wx(n),wy(n)

      wx=x
      wy=y
      do i=1,n
        in=n-i+1
        x(i)=wx(in)
        y(i)=wy(in)
      enddo

      return
      end
