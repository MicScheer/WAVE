*CMZ :          11/04/2018  12.06.39  by  Michael Scheer
*-- Author :    Michael Scheer   06/04/2018
      subroutine util_rotate_matrix(axis,angle,rotmat,istat)

      implicit none

      double precision axis(3),angle,rotmat(3,3),
     &  ex(3),ey(3),ez(3),
     &  erx(3),ery(3),erz(3)

      integer istat,kstat

      data ex/1.0d0,0.0d0,0.0d0/
      data ey/0.0d0,1.0d0,0.0d0/
      data ez/0.0d0,0.0d0,1.0d0/

      istat=0

      if (angle.eq.0.0d0) then
        rotmat=0.0d0
        rotmat(1,1)=1.0d0
        rotmat(2,2)=1.0d0
        rotmat(3,3)=1.0d0
        return
      endif

      call util_rotate_vector(ex,axis,angle,erx,kstat)
      if (kstat.ne.0) istat=-1
      call util_rotate_vector(ey,axis,angle,ery,kstat)
      if (kstat.ne.0) istat=-1
      call util_rotate_vector(ez,axis,angle,erz,kstat)
      if (kstat.ne.0) istat=-1

      rotmat(1:3,1)=erx(1:3)
      rotmat(1:3,2)=ery(1:3)
      rotmat(1:3,3)=erz(1:3)

      return
      end
