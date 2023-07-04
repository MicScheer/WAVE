*CMZ :  3.05/28 07/01/2019  13.00.10  by  Michael Scheer
*CMZ :  3.05/26 07/12/2018  12.07.54  by  Michael Scheer
*CMZ :  3.05/20 01/11/2018  15.08.03  by  Michael Scheer
*CMZ : 00.00/19 07/06/2016  12.17.28  by  Michael Scheer
*-- Author :    Michael Scheer   07/06/2016
      subroutine util_plane_hit_hesse(pp,vp,pl,vl,hit,dist,distn,istat)

      implicit none

      double precision pp(3),vp(3),pl(3),vl(3),hit(3),p1(3),p2(3),
     &  u(3),v(3),w(3),vln,vpn,dist,distn

      integer istat

      vpn=sqrt(vp(1)**2+vp(2)**2+vp(3)**2)
      if (vpn.eq.0.0d0) then
        istat=-1
        return
      endif

      vln=sqrt(vl(1)**2+vl(2)**2+vl(3)**2)
      if (vln.eq.0.0d0) then
        istat=-2
        return
      endif

      u=vp/vpn

      if (abs(u(3)).lt.0.5d0) then
        v(1)=-u(2)
        v(2)=u(1)
        v(3)=0.0d0
      else
        v(1)=0.0d0
        v(2)=-u(3)
        v(3)=u(2)
      endif

      w(1)=u(2)*v(3)-u(3)*v(2)
      w(2)=u(3)*v(1)-u(1)*v(3)
      w(3)=u(1)*v(2)-u(2)*v(1)

      p1=pp+v
      p2=pp+w

      call util_plane_hit(pp,p1,p2,pl,vl,hit,dist,distn,u,istat)

      return
      end
