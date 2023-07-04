*CMZ :          08/08/2020  17.24.36  by  Michael Scheer
*-- Author :    Michael Scheer   08/08/2020
      subroutine util_check_smoothness(n,x,y,ys,threshhold,is,mode,istat)
      implicit none

      double precision x(n),y(n),ys(n),threshhold,ymin,ymax,thresh,g
      double  precision, dimension (:), allocatable :: xa,ya,y2a
      double precision :: tiny=1.0e-12

      integer n,is(n),istat,mode,i,kbad,k,i1,i2,lgood

      ! check smoothness and fix guessed bad points
      ! first and last points are assumed to be good

      !mode:
      !0: threshhold is relative to dy
      !1: threshhold is relative to max(ymax-ymin)
      !else: threshhold is absolute

      !output:
      !ys: smoothed y(x)
      !is: is(i)=-1: assumed to be bad
      !is: is(i)=+1: assumed to be fixed
      !    is(i)=0: assumed to be good
      !istat: Error if negative or number of bad points

      is=0
      kbad=0

      if (mode.eq.1) then
        ymin=1.0d30
        ymax=-1.0d30
        do i=1,n
          if (y(i).lt.ymin) ymin=y(i)
          if (y(i).gt.ymax) ymax=y(i)
        enddo
        thresh=(ymax-ymin)*abs(threshhold)
      else if (mode.ne.0.and.mode.ne.1) then
        thresh=abs(threshhold)
      endif

      istat=0
      lgood=1

      do i=2,n-1
        i1=lgood
        i2=i+1
        g=y(i1)+(y(i2)-y(i1))/(x(i2)-x(i1))*(x(i)-x(i1))
        if (mode.eq.0) then
          thresh=abs(y(i2)-y(i1))*threshhold
        endif
        if (abs(g-y(i)).gt.thresh) then
          is(i)=-1
          kbad=1
        else
          lgood=i
        endif
      enddo

      if (kbad.ne.0) then

        k=0

        if (is(1).ne.0.and.is(2).ne.0)
     &    kbad=kbad+1
        if (is(n-1).ne.0.and.is(n).ne.0)
     &    kbad=kbad+1

        allocate(xa(n),ya(n),y2a(n))

        if (is(1).eq.0) then
          k=k+1
          xa(k)=x(1)
          ya(k)=y(1)
        endif

        do i=2,n-1
          if (is(i).eq.0) then
            k=k+1
            xa(k)=x(i)
            ya(k)=y(i)
          endif
          if (is(i).ne.0.and.(is(i-1).ne.0.or.is(i+1).ne.0)) kbad=kbad+1
        enddo

        if (is(n).eq.0) then
          k=k+1
          xa(k)=x(n)
          ya(k)=y(n)
        endif

        do i=1,n
          if (is(i).eq.-1) then
            call util_spline_interpolation_f90(k,xa,ya,x(i),ys(i),y2a,0)
          else
            ys(i)=y(i)
          endif
        enddo

        deallocate(xa,ya,y2a)

        istat=0

        do i=1,n
          if (is(i).eq.-1) then
            if (mode.eq.0) then
              thresh=abs(ys(i2)-ys(i1))*threshhold
            endif
            if (ys(i)+y(i).ne.0.0d0
     &          .and.abs((ys(i)-y(i))/(ys(i)+y(i))).lt.tiny) then
              is(i)=0
            endif
          endif !(is(i).eq.-1) then
          if (is(i).lt.0) istat=istat+1
        enddo

      else

        ys=y

      endif

      return
      end
