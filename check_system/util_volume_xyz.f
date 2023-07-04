*CMZ :          06/03/2023  21.48.22  by  Michael Scheer
*CMZ :  2.04/03 03/03/2023  15.00.22  by  Michael Scheer
*CMZ :  2.04/02 27/02/2023  18.38.39  by  Michael Scheer
*-- Author :    Michael Scheer   26/02/2023
      subroutine util_volume_xyz(n,xin,yin,zin,tiny,v,kfail)

      implicit none

      double precision :: tiny

      double precision xin(*),yin(*),zin(*),v,a,vnor(3),dist,gc(3)

      integer l,n,kfail,nhull,nedge,nface,ipoi,iface,npoi,kfacelast,k,iover,i

      double precision, dimension(:), allocatable :: x,y,z,hull
      integer, dimension(:,:), allocatable :: kedge
      integer, dimension(:), allocatable :: kface,khull

      ! khull(n)
      ! kedge(4,2*n-2)
      ! kface((n+1)*n)

      allocate(x(n),y(n),z(n),hull(n),khull(n),kedge(4,2*n-2),kface((n+1)*n))

      call util_convex_hull_3d(
     &  n,xin,yin,zin,khull,kedge,kface,nhull,nedge,nface,kfacelast,tiny,kfail)

      v=0.0d0
      if (kfail.ne.0.or.nhull.le.3) goto 9999

      gc=0.0d0

      do i=1,nhull
        k=khull(i)
        gc=gc+[xin(k),yin(k),zin(k)]
      enddo
      gc=gc/dble(n)

      l=0
      do iface=1,nface
        npoi=kface(l+1)
        do ipoi=1,npoi
          k=kface(l+1+ipoi)
          x(ipoi)=xin(k)
          y(ipoi)=yin(k)
          z(ipoi)=zin(k)
        enddo
        call util_plane(
     &    [x(1),y(1),z(1)],[x(2),y(2),z(2)],[x(3),y(3),z(3)],
     &    gc,vnor,dist,iover,kfail)
        if (kfail.ne.0) then
          v=0.0d0
          goto 9999
        endif
        call util_area_xyz(npoi,x,y,z,tiny,a,kfail)
        if (kfail.ne.0) then
          v=0.0d0
          goto 9999
        endif
        v=v+a*abs(dist)
        l=l+npoi+1
      enddo

      v=v/3.0d0

9999  deallocate(x,y,z,hull,kedge,kface,khull)

      return
      end
