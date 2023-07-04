*CMZ :  3.05/23 28/11/2018  09.28.33  by  Michael Scheer
*-- Author :    Michael Scheer   08/03/2018
      subroutine util_rotate_vector_to_x_axis(vin,rotmat,istat)

      implicit none

      double precision cosphi,sinphi,costhe,sinthe,vin(3),
     &  vx,vy,vz,vn,rotmat(3,3),rotphi(3,3),rotthe(3,3)

      integer istat

      istat=0
      vn=sqrt(vin(1)**2+vin(2)**2+vin(3)**2)

      rotmat=0.0d0

      if (vn.eq.0.0d0) then
        istat=-1
        return
      endif

      vx=vin(1)/vn
      vy=vin(2)/vn
      vz=vin(3)/vn

      if (vy**2+vz**2.lt.1.0d-9) then
        rotmat(1,1)=1.0d0
        rotmat(1,2)=0.0d0
        rotmat(1,3)=0.0d0
        rotmat(2,1)=0.0d0
        rotmat(2,2)=1.0d0
        rotmat(2,3)=0.0d0
        rotmat(3,1)=0.0d0
        rotmat(3,2)=0.0d0
        rotmat(3,3)=1.0d0
        return
      endif

      costhe=vz
      sinthe=sqrt(1.0d0-min(1.0d0,costhe**2))

      if (sinthe.lt.1.0d-9) then
        rotmat(1,1)=0.0d0
        rotmat(1,2)=0.0d0
        rotmat(1,3)=1.0d0
        rotmat(2,1)=0.0d0
        rotmat(2,2)=1.0d0
        rotmat(2,3)=0.0d0
        rotmat(3,1)=-1.0d0
        rotmat(3,2)=0.0d0
        rotmat(3,3)=0.0d0
        return
      endif

      if (vz.lt.0.0d0) sinthe=-sinthe

      if (abs(vz).lt.1.0d-9) then
        cosphi=vx
        sinphi=vy
        rotmat(1,1)=cosphi
        rotmat(1,2)=sinphi
        rotmat(1,3)=0.0d0
        rotmat(2,1)=-sinphi
        rotmat(2,2)=cosphi
        rotmat(2,3)=0.0d0
        rotmat(3,1)=0.0d0
        rotmat(3,2)=0.0d0
        rotmat(3,3)=1.0d0
        return
      endif

      cosphi=vx/sinthe
      sinphi=sqrt(1.0d0-min(1.0d0,cosphi**2))
      if(vy.lt.0.0d0) sinphi=-sinphi

      rotphi(1,1)=cosphi
      rotphi(1,2)=sinphi
      rotphi(1,3)=0.0d0
      rotphi(2,1)=-sinphi
      rotphi(2,2)=cosphi
      rotphi(2,3)=0.0d0
      rotphi(3,1)=0.0d0
      rotphi(3,2)=0.0d0
      rotphi(3,3)=1.0d0

      rotthe(1,1)=sinthe
      rotthe(1,2)=0.0d0
      rotthe(1,3)=costhe

      rotthe(2,1)=0.0
      rotthe(2,2)=1.0d0
      rotthe(2,3)=0.0

      rotthe(3,1)=-costhe
      rotthe(3,2)=0.0d0
      rotthe(3,3)=sinthe

      rotmat(1,1)=
     &  rotthe(1,1)*rotphi(1,1)+rotthe(1,2)*rotphi(2,1)+rotthe(1,3)*rotphi(3,1)
      rotmat(1,2)=
     &  rotthe(1,1)*rotphi(1,2)+rotthe(1,2)*rotphi(2,2)+rotthe(1,3)*rotphi(3,2)
      rotmat(1,3)=
     &  rotthe(1,1)*rotphi(1,3)+rotthe(1,2)*rotphi(2,3)+rotthe(1,3)*rotphi(3,3)

      rotmat(2,1)=
     &  rotthe(2,1)*rotphi(1,1)+rotthe(2,2)*rotphi(2,1)+rotthe(2,3)*rotphi(3,1)
      rotmat(2,2)=
     &  rotthe(2,1)*rotphi(1,2)+rotthe(2,2)*rotphi(2,2)+rotthe(2,3)*rotphi(3,2)
      rotmat(2,3)=
     &  rotthe(2,1)*rotphi(1,3)+rotthe(2,2)*rotphi(2,3)+rotthe(2,3)*rotphi(3,3)

      rotmat(3,1)=
     &  rotthe(3,1)*rotphi(1,1)+rotthe(3,2)*rotphi(2,1)+rotthe(3,3)*rotphi(3,1)
      rotmat(3,2)=
     &  rotthe(3,1)*rotphi(1,2)+rotthe(3,2)*rotphi(2,2)+rotthe(3,3)*rotphi(3,2)
      rotmat(3,3)=
     &  rotthe(3,1)*rotphi(1,3)+rotthe(3,2)*rotphi(2,3)+rotthe(3,3)*rotphi(3,3)

      return
      end
