*CMZ :          30/05/2019  21.20.00  by  Michael Scheer
*-- Author :    Michael Scheer   12/04/2016
      subroutine util_weed_triangle(n,x,y,itri,l,irest)

      implicit none

      double precision x(*),y(*),xx1,yy1,xx2,yy2,xx3,yy3,
     &  x1,y1,x2,y2,x3,y3,dx1,dy1,dx2,dy2,dx3,dy3

      integer n,itri(3),l(0:n),irest(0:n),i

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
      dx3=x1-x3
      dy3=y1-y3

      irest(0)=0
      do i=1,l(0)
        if(l(i).eq.itri(1).or.l(i).eq.itri(2).or.l(i).eq.itri(3)) cycle
        xx1=x(l(i))-x1
        yy1=y(l(i))-y1
        xx2=x(l(i))-x2
        yy2=y(l(i))-y2
        xx3=x(l(i))-x3
        yy3=y(l(i))-y3
        if (
     &      (xx1*dy1-yy1*dx1.lt.0.0d0.or.
     &      xx2*dy2-yy2*dx2.lt.0.0d0.or.
     &      xx3*dy3-yy3*dx3.lt.0.0d0)
     &      .and.
     &      (xx1*dy1-yy1*dx1.gt.0.0d0.or.
     &      xx2*dy2-yy2*dx2.gt.0.0d0.or.
     &      xx3*dy3-yy3*dx3.gt.0.0d0)
     &    ) then
          irest(0)=irest(0)+1
          irest(irest(0))=l(i)
        endif
      enddo

      return
      end
