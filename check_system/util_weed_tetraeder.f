*CMZ :          07/06/2019  14.43.27  by  Michael Scheer
*-- Author :    Michael Scheer   04/06/2019
      subroutine util_weed_tetraeder(n,xb,yb,zb,ibuff,itet,irest,tiny)

      implicit none

      double precision xb(n),yb(n),zb(n),dist,sp,tiny,
     &  dpg(3),p1(3),p2(3),p3(3),p4(3),q(3),vnor(3),g(3)

      integer i,k,n,itet(4),irest(0:n),in,ibuff(0:n),iover,istat

      p1(1)=xb(itet(1))
      p1(2)=yb(itet(1))
      p1(3)=zb(itet(1))

      p2(1)=xb(itet(2))
      p2(2)=yb(itet(2))
      p2(3)=zb(itet(2))

      p3(1)=xb(itet(3))
      p3(2)=yb(itet(3))
      p3(3)=zb(itet(3))

      p4(1)=xb(itet(4))
      p4(2)=yb(itet(4))
      p4(3)=zb(itet(4))

      g=(p1+p2+p3+p4)/4.0d0

      irest=0
      do k=1,ibuff(0)

        i=ibuff(k)

        if (i.eq.itet(1).or.i.eq.itet(2).or.i.eq.itet(3).or.i.eq.itet(4))
     &    cycle

        in=0

        q(1)=xb(i)
        q(2)=yb(i)
        q(3)=zb(i)

        dpg=p1-g
        call util_plane(p1,p2,p3,q,vnor,dist,iover,istat)
        sp=dpg(1)*vnor(1)+dpg(2)*vnor(2)+dpg(3)*vnor(3)
        if (abs(dist).le.tiny.or.
     &    sp.gt.0.0d0.and.dist.le.0.0d0.or.
     &    sp.lt.0.0d0.and.dist.ge.0.0d0
     &    ) in=in+1

        dpg=p1-g
        call util_plane(p1,p2,p4,q,vnor,dist,iover,istat)
        sp=dpg(1)*vnor(1)+dpg(2)*vnor(2)+dpg(3)*vnor(3)
        if (abs(dist).le.tiny.or.
     &    sp.gt.0.0d0.and.dist.le.0.0d0.or.
     &    sp.lt.0.0d0.and.dist.ge.0.0d0
     &    ) in=in+1

        dpg=p2-g
        call util_plane(p2,p3,p4,q,vnor,dist,iover,istat)
        sp=dpg(1)*vnor(1)+dpg(2)*vnor(2)+dpg(3)*vnor(3)
        if (abs(dist).le.tiny.or.
     &    sp.gt.0.0d0.and.dist.le.0.0d0.or.
     &    sp.lt.0.0d0.and.dist.ge.0.0d0
     &    ) in=in+1

        dpg=p3-g
        call util_plane(p3,p1,p4,q,vnor,dist,iover,istat)
        sp=dpg(1)*vnor(1)+dpg(2)*vnor(2)+dpg(3)*vnor(3)
        if (abs(dist).le.tiny.or.
     &    sp.gt.0.0d0.and.dist.le.0.0d0.or.
     &    sp.lt.0.0d0.and.dist.ge.0.0d0
     &    ) in=in+1

        if (in.eq.4) cycle

        irest(0)=irest(0)+1
        irest(irest(0))=i

      enddo

      return
      end
