*CMZ :          06/06/2019  16.34.47  by  Michael Scheer
*CMZ : 00.00/19 07/06/2016  12.17.28  by  Michael Scheer
*-- Author :    Michael Scheer   07/06/2016
      subroutine util_plane_npoints(n,poi,p,vnor,dist,g,iover,istat)

      implicit none

      double precision poi(3,n),p1(3),p2(3),p3(3),p(3),vnor(3),dist,
     &  v21(3),vin(3),v32(3),vn,g(3),vpg(3),vc(3),dp(3)

      integer iover,istat,n,i

      istat=0

      g=0.0d0
      do i=1,n
        g=g+poi(1:3,i)
      enddo
      g=g/dble(n)


      p1=poi(1:3,1)
      p2=poi(1:3,2)
      p3=poi(1:3,3)

      v21=p2-p1
      v32=p3-p2

      vnor(1)=v21(2)*v32(3)-v21(3)*v32(2)
      vnor(2)=v21(3)*v32(1)-v21(1)*v32(3)
      vnor(3)=v21(1)*v32(2)-v21(2)*v32(1)
      vn=sqrt(vnor(1)**2+vnor(2)**2+vnor(3)**2)

      if (vn.le.0.0d0) then
        istat=-1
        return
      endif

      vnor=vnor/vn

      vpg=p-g
      dist=vpg(1)*vnor(1)+vpg(2)*vnor(2)+vpg(3)*vnor(3)

      iover=1

      do i=1,n
        if (i.lt.n) then
          dp=poi(:,i+1)-poi(:,i)
        else
          dp=poi(:,1)-poi(:,n)
        endif
        vin=p-poi(:,i)-dist*vnor
        vc(1)=dp(2)*vin(3)-dp(3)*vin(2)
        vc(2)=dp(3)*vin(1)-dp(1)*vin(3)
        vc(3)=dp(1)*vin(2)-dp(2)*vin(1)
        if (vc(1)*vnor(1)+vc(2)*vnor(2)+vc(3)*vnor(3).lt.0.0d0) then
          iover=0
          exit
        endif
      enddo

      return
      end
