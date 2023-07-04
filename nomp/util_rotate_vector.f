*CMZ :  3.05/23 11/04/2018  12.04.00  by  Michael Scheer
*-- Author :    Michael Scheer   06/04/2018
      subroutine util_rotate_vector(vin,axis,angle,vout,istat)

      implicit none

      double precision vin(3),vout(3),axis(3),angle,vn(3),vnor,
     &  anor,ex(3),ey(3),ez(3),
     &  cosz,cosx

      integer istat

      istat=0

      vnor=sqrt(vin(1)**2+vin(2)**2+vin(3)**2)

      if (vnor.eq.0.0d0) then
        istat=1
        vout=vin
        return
      endif

      anor=sqrt(axis(1)**2+axis(2)**2+axis(3)**2)

      if (anor.eq.0.0d0) then
        istat=2
        vout=vin
        return
      endif

      if (angle.eq.0.0d0) then
        vout=vin
        return
      endif

      vn=vin/vnor
      ez=axis/anor

      ey(1)=ez(2)*vn(3)-ez(3)*vn(2)
      ey(2)=ez(3)*vn(1)-ez(1)*vn(3)
      ey(3)=ez(1)*vn(2)-ez(2)*vn(1)

      if (ey(1).eq.0.0d0.and.ey(2).eq.0.0d0.and.ey(3).eq.0.0d0) then
        vout=vin
        return
      endif

      ey=ey/sqrt(ey(1)**2+ey(2)**2+ey(3)**2)

      ex(1)=ey(2)*ez(3)-ey(3)*ez(2)
      ex(2)=ey(3)*ez(1)-ey(1)*ez(3)
      ex(3)=ey(1)*ez(2)-ey(2)*ez(1)

      cosx=vin(1)*ex(1)+vin(2)*ex(2)+vin(3)*ex(3)
      cosz=vin(1)*ez(1)+vin(2)*ez(2)+vin(3)*ez(3)

      vout=cosx*(cos(angle)*ex+sin(angle)*ey)+cosz*ez

      return
      end
