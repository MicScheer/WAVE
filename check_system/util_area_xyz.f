*CMZ :          07/03/2023  09.33.46  by  Michael Scheer
*-- Author :    Michael Scheer   25/02/2023
      subroutine util_area_xyz(n,x,y,z,tiny,a,kfail)
      implicit none

      double precision :: x(*),y(*),z(*),a,tiny,p1(3),p2(3),p3(3),vnor(3),vn
      integer n,i,kfail,nh,k,i1

      integer, dimension (:), allocatable :: ihull

      a=0.0d0

      if (n.lt.3) then
        kfail=-1
        return
      endif

      allocate(ihull(n+1))

      call util_convex_hull_2d(n,x,y,nh,ihull,tiny,kfail)

      if (nh.lt.4) kfail=1

      if (kfail.eq.0) then

        i=ihull(1)
        p1=[x(i),y(i),z(i)]
        i=ihull(2)
        p2=[x(i),y(i),z(i)]
        i=ihull(3)
        p3=[x(i),y(i),z(i)]

        call util_vcross(p2-p1,p3-p2,vnor)
        vn=sqrt(vnor(1)**2+vnor(2)**2+vnor(3)**2)

        vnor=abs(vnor)/vn

        if (vnor(3).gt.0.1d0) then
          do k=1,nh-1
            i=ihull(k)
            i1=ihull(k+1)
            a=a+(y(i)+y(i1))/2.0d0*(x(i1)-x(i))/vnor(3)
          enddo
        else
          kfail=1
        endif

      endif

      if (kfail.ne.0) then

        call util_convex_hull_2d(n,y,z,nh,ihull,tiny,kfail)
        if (nh.lt.4) kfail=1

        if (kfail.eq.0) then

          i=ihull(1)
          p1=[x(i),y(i),z(i)]
          i=ihull(2)
          p2=[x(i),y(i),z(i)]
          i=ihull(3)
          p3=[x(i),y(i),z(i)]

          call util_vcross(p2-p1,p3-p2,vnor)
          vn=sqrt(vnor(1)**2+vnor(2)**2+vnor(3)**2)

          vnor=abs(vnor)/vn

          if (vnor(1).gt.0.1d0) then
            do k=1,nh-1
              i=ihull(k)
              i1=ihull(k+1)
              a=a+(z(i)+z(i1))/2.0d0*(y(i1)-y(i))/vnor(1)
            enddo
          else
            kfail=1
          endif

        endif
      endif

      if (kfail.ne.0) then

        call util_convex_hull_2d(n,x,z,nh,ihull,tiny,kfail)
        if (nh.lt.4) kfail=1

        if (kfail.eq.0) then

          i=ihull(1)
          p1=[x(i),y(i),z(i)]
          i=ihull(2)
          p2=[x(i),y(i),z(i)]
          i=ihull(3)
          p3=[x(i),y(i),z(i)]

          call util_vcross(p2-p1,p3-p2,vnor)
          vn=sqrt(vnor(1)**2+vnor(2)**2+vnor(3)**2)

          vnor=abs(vnor)/vn

          if (vnor(2).gt.0.1d0) then
            do k=1,nh-1
              i=ihull(k)
              i1=ihull(k+1)
              a=a+(z(i)+z(i1))/2.0d0*(x(i1)-x(i))/vnor(2)
            enddo
          else
            kfail=1
          endif

        endif
      endif

      deallocate(ihull)

      a=abs(a)

      return
      end
