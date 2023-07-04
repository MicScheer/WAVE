*CMZ : 00.00/14 17/09/2011  20.11.16  by  Michael Scheer
*CMZ : 00.00/11 09/03/2011  15.33.07  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_pspline_minmax(n,xa,ya,xmin,ymin,xmax,ymax,ifail)

      implicit none

      double precision, dimension (:), allocatable :: y1

      double precision xa(n),ya(n),yp,ypp,ymin,xmin,ymax,xmax,
     &  a(3),xzero(2),y1p(3),xopt,yopt,yzero

      integer n,ifail,ifound,i,i1

      allocate (y1(n))

      ymin=1.0d30
      ymax=-1.0d30

      if (n.lt.3) then
        ifail=1
        return
      endif

      do i=1,n
        if (ya(i).lt.ymin) then
          xmin=xa(i)
          ymin=ya(i)
        endif
        if (ya(i).gt.ymax) then
          xmax=xa(i)
          ymax=ya(i)
        endif
      enddo

      call util_spline_coef_second_order(xa,ya,n,9999.0d0,1,y1,ifail)
      if (ifail.ne.0) return

      ifail=-2
      ifound=0

      do i=1,n-1

        if (
     &      y1(i).le.0.0d0.and.y1(i+1).gt.0.0d0
     &      .or.
     &      y1(i).gt.0.0d0.and.y1(i+1).le.0.0d0) then
          ifound=ifound+1
          if(i.gt.1.and.i.lt.n) then
            i1=i-1
          else if(i.eq.1) then
            i1=i
          else if(i.eq.n) then
            i1=n-2
          endif
          call util_parabel_zero(xa(i1),y1(i1),a,y1p,xopt,yopt,xzero,ifail)

          if (ifail.ne.0) return

          ifail=-2

          if (
     &        xa(i1).lt.xa(i1+2).and.xzero(1).ge.xa(i1).and.xzero(1).le.xa(i1+2)
     &        .or.
     &        xa(i1).gt.xa(i1+2).and.xzero(1).le.xa(i1).and.xzero(1).ge.xa(i1+2)
     &        ) then
            call util_spline_inter_second_order(xa,ya,y1,n,xzero(1),yzero,
     &        yp,ypp,-1,ifail)
            if (ifail.ne.0) return
            if(yzero.gt.ymax) then
              xmax=xzero(1)
              ymax=yzero
            endif
            if(yzero.lt.ymin) then
              xmin=xzero(1)
              ymin=yzero
            endif
            ifail=0
          endif

          if (
     &        xa(i1).lt.xa(i1+2).and.xzero(2).ge.xa(i1).and.xzero(2).le.xa(i1+2)
     &        .or.
     &        xa(i1).gt.xa(i1+2).and.xzero(2).le.xa(i1).and.xzero(2).ge.xa(i1+2)
     &        ) then
            call util_spline_inter_second_order(xa,ya,y1,n,xzero(2),yzero,
     &        yp,ypp,-1,ifail)
            if (ifail.ne.0) return
            if(yzero.gt.ymax) then
              xmax=xzero(2)
              ymax=yzero
            endif
            if(yzero.lt.ymin) then
              xmin=xzero(2)
              ymin=yzero
            endif
            ifail=0
          endif
        endif

      enddo

      deallocate (y1)

      return
      end
