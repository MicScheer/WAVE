*CMZ :          25/02/2023  16.23.46  by  Michael Scheer
*-- Author :    Michael Scheer   25/02/2023
      subroutine util_convex_hull_3d_overwrite(
     &  nin,xin,yin,zin,khull,kedge,kface,
     &  nhull,nedge,nface,kfacelast,tiny,
     &  kfail)

      ! like  util_convex_hull_3d, but overwrites xin,yin,zin

      implicit none

      double precision xin(*),yin(*),zin(*),tiny

      double precision, dimension (:), allocatable ::  xh,yh,zh

      integer khull(*),kedge(4,*),kface(*)
      integer nin,n,k,kfail,i,nhull,nedge,nface,kfacelast

      call util_convex_hull_3d(
     &  nin,xin,yin,zin,khull,kedge,kface,
     &  nhull,nedge,nface,kfacelast,tiny,
     &  kfail)

      if (kfail.ne.0) then
        return
      endif

      allocate(xh(nhull),yh(nhull),zh(nhull))

      n=nhull
      do i=1,n
        k=khull(i)
        xh(i)=xin(k)
        yh(i)=yin(k)
        zh(i)=zin(k)
      enddo

      call util_convex_hull_3d(
     &  n,xh,yh,zh,khull,kedge,kface,
     &  nhull,nedge,nface,kfacelast,tiny,
     &  kfail)

      n=nhull
      do i=1,n
        xin(i)=xh(i)
        yin(i)=yh(i)
        zin(i)=zh(i)
      enddo

      deallocate(xh,yh,zh)

      return
      end
