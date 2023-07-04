*CMZ :          07/06/2019  15.48.29  by  Michael Scheer
*-- Author :    Michael Scheer   04/06/2019
      subroutine util_next_point_of_hull(n,xb,yb,zb,itri,irest,next)

      implicit none

      double precision xb(n),yb(n),zb(n),
     &  p1(3),p2(3),p3(3),p(3),vnor(3),dist,distmax

      integer i,n,k,itri(3),irest(0:n),next,iover,istat

      p1(1)=xb(itri(1))
      p1(2)=yb(itri(1))
      p1(3)=zb(itri(1))

      p2(1)=xb(itri(2))
      p2(2)=yb(itri(2))
      p2(3)=zb(itri(2))

      p3(1)=xb(itri(3))
      p3(2)=yb(itri(3))
      p3(3)=zb(itri(3))

      distmax=0.0d0
      next=0

      do i=1,irest(0)

        k=irest(i)

        p(1)=xb(k)
        p(2)=yb(k)
        p(3)=zb(k)

        call util_plane(p1,p2,p3,p,vnor,dist,iover,istat)

        if (iover.eq.1.and.dist.gt.distmax) then
          next=k
          distmax=dist
        endif

      enddo

      return
      end
