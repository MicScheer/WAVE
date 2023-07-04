*CMZ : 00.00/11 07/06/2011  14.38.25  by  Michael Scheer
*-- Author :    Michael Scheer   29/03/2011
      subroutine util_average_two_functions(n1,x1,y1,n2,x2,y2,
     &  weight1,weight2,x,y,cmode)

      implicit none

      double precision, dimension (:), allocatable :: ws1,ws2,ws3,ws4,coef

      double precision x1(n1),y1(n1),x2(n2),y2(n2),weight1,weight2,x(*),y(*),
     &  dist,dist0,distmin,w12,yy
      character(2)cmode

      integer n1,n2,i1,i2,last

      w12=1.0d0
      if (abs(weight1)+abs(weight2).gt.0.0d0)
     &  w12=1.0d0/(abs(weight1)+abs(weight2))

      if (cmode.eq.'hv') then

        last=1
        dist0=1.0d30

        do i1=1,n1
          distmin=1.0d30
          dist0=1.0d30
          do i2=last,n2
            dist=(x1(i1)-x2(i2))**2+(y1(i1)-y2(i2))**2
            if (dist.lt.distmin) then
              distmin=dist
              last=i2
            endif
            if (dist.ge.dist0) goto 9
            dist0=dist
          enddo
9         continue

          x(i1)=(weight1*x1(i1)+weight2*x2(last))*w12
          y(i1)=(weight1*y1(i1)+weight2*y2(last))*w12

        enddo

      else !cmode

        allocate(ws1(n2))
        allocate(ws2(n2))
        allocate(ws3(n2))
        allocate(ws4(n2))
        allocate(coef(n2))

        call util_spline_coef(x2,y2,n2,9999.0d0,9999.0d0,coef,ws1,ws2,ws3,ws4)

        x(1)=x1(1)
        call util_spline_inter(x2,y2,coef,n2,x1(1),yy,-1)
        y(1)=(weight1*y1(1)+weight2*yy)*w12

        do i1=2,n1
          x(i1)=x1(i1)
          call util_spline_inter(x2,y2,coef,n2,x1(i1),yy,0)
          y(i1)=(weight1*y1(i1)+weight2*yy)*w12
        enddo

        deallocate(ws1)
        deallocate(ws2)
        deallocate(ws3)
        deallocate(ws4)
        deallocate(coef)

      endif !cmode

      return
      end
