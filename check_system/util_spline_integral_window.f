*CMZ : 00.00/16 18/03/2014  17.02.27  by  Michael Scheer
*CMZ : 00.00/15 12/10/2013  12.18.16  by  Michael Scheer
*CMZ : 00.00/07 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/02 17/08/2004  09.47.26  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_spline_integral_window(x,y,n,xmin,xmax,resultat
     &  ,coef,work1,work2,work3,work4,mode,istat)

c---  calculates intergral of y(x) in interval [xmin,xmax] via splines

      implicit none

      real*8 x(n),y(n),resultat,xmin,xmax,ymin,ymax,cmin,cmax
      real*8 coef(n),work1(n),work2(n),work3(n),work4(n)

      integer i,n,mode,imn,imx,i1,i2,istat

      resultat=0.0d0
      istat=-1

      if (x(1).gt.x(n)) return

c---  spline-coefficients

      if (mode.lt.0) then
        call util_spline_coef(x,y,n,9999.0d0,9999.0d0,coef,work1,work2,work3,work4)
      endif

      if (xmin.lt.x(1)) then
        istat=-2
        return
      endif
      if (xmax.gt.x(n)) then
        istat=-3
        return
      endif

c intervals

      call util_bisec(n,x,xmin,imn)
      call util_bisec(n,x,xmax,imx)

      if (imn.eq.imx) then

c from xmin to xmax

        call util_spline_inter(x,y,coef,n,xmin,ymin,0)
        call util_spline_inter(x,y,coef,n,xmax,ymax,0)

        cmin=coef(imn)+(coef(imn+1)-coef(imn))/(x(imn+1)-x(imn))*(xmin-x(imn))
        cmax=coef(imn)+(coef(imn+1)-coef(imn))/(x(imn+1)-x(imn))*(xmax-x(imn))

        resultat=resultat
     &    +(xmax-xmin)*0.5d0
     &    *(ymin+ymax)
     &    -(xmax-xmin)**3/24.d0
     &    *(cmin+cmax)

        return

      endif

      i1=imn+1
      i2=imx-1

c--- integration

c from xmin to x(imn+1)
      call util_spline_inter(x,y,coef,n,xmin,ymin,0)
      cmin=coef(imn)+(coef(imn+1)-coef(imn))/(x(imn+1)-x(imn))*(xmin-x(imn))
      resultat=resultat
     &  +(x(imn+1)-xmin)*0.5d0
     &  *(ymin+y(imn+1))
     &  -(x(imn+1)-xmin)**3/24.d0
     &  *(cmin+coef(imn+1))

c from x(imx) to xmax
      call util_spline_inter(x,y,coef,n,xmax,ymax,0)
      cmax=coef(imx)+(coef(imx+1)-coef(imx))/(x(imx+1)-x(imx))*(xmax-x(imx))
      resultat=resultat+
     &  (xmax-x(imx))*0.5d0
     &  *(ymax+y(imx))
     &  -(xmax-x(imx))**3/24.d0
     &  *(cmax+coef(imx))

      do i=i1,i2

        resultat=resultat
     &    +(x(i+1)-x(i))*0.5d0
     &    *(y(i)+y(i+1))
     &    -(x(i+1)-x(i))**3/24.d0
     &    *(coef(i)+coef(i+1))

      enddo

      istat=0

      return
      end
