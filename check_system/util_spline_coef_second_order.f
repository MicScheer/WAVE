*CMZ : 00.00/09 17/12/2010  09.45.53  by  Michael Scheer
*CMZ : 00.00/07 12/10/2009  12.17.45  by  Michael Scheer
*CMZ : 00.00/02 14/04/2003  12.46.09  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.48  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE util_spline_coef_second_order(X,Y,N,YP1,nyp,Y1,ifail)

C--- CALCULATES SPLINE COEFFICIENTS for a spline up to second order

C--   INPUT:

C-       N: NUMBER OF X,Y-VALUES
C-       X: ARRAY OF X-VALUES
C-       Y: ARRAY OF Y-VALUES
C-       YP1:  first DERIVATIVE AT FIRST NYPth-VALUE

C--   OUPUT:

C-       Y1:   SPLINE-COEFFICIENTS
c-    ifail: Error status, should be zero

c fl(x)=yl+yp1l*(x-xl)+1/2*yppl*(x-xl)**2
c flp(x)=yp1l+yppl*(x-xl)
c fl(xh)=yl+yp1l*(xh-xl)+1/2*yppl*(xh-xl)**2 = yh
c flp(x)=yp1l+yppl*(x-xl) = yp1h
c
c yppl = 2 * (yh - yl -yp1l*(xh-xl)) / (xh-xl)**2

      IMPLICIT NONE

      INTEGER N,J,nyp
      REAL*8  X(N),Y(N),Y1(N),yp1,yppl

      double precision xx(3),yy(3),a(3),yp(3),xopt,yopt
      INTEGER ifail

      ifail=-1

      if (nyp.lt.0.or.nyp.gt.n) return
      if (n.le.1) return

      if (n.eq.2) then

        if(x(2).eq.x(1)) return

        if(yp1.eq.9999.0d0) then
          y1(1)=(y(2)-y(1))/(x(2)-x(1))
          y1(2)=y1(1)
          ifail=0
          return
        endif

        if (nyp.eq.1) then
          y1(1)=yp1
          yppl=2.0d0*(y(2)-y(1)-y1(1)*(x(2)-x(1)))/(x(2)-x(1))**2
          y1(2)=y1(1)+yppl*(x(2)-x(1))
        else if (nyp.eq.2) then
          y1(2)=yp1
          yppl=2.0d0*(y(1)-y(2)-y1(2)*(x(1)-x(2)))/(x(1)-x(2))**2 ! exchanged
          y1(1)=y1(2)+yppl*(x(1)-x(2))
        endif

        ifail=0
        return
      endif

      y1(nyp)=yp1

      if (abs(yp1).eq.9999.0d0) then

        if (nyp.eq.1) then
          xx(1:3)=x(1:3)
          yy(1:3)=y(1:3)
        else if (nyp.eq.n) then
          xx(1:3)=x(n-2:n)
          yy(1:3)=y(n-2:n)
        else
          xx(1:3)=x(nyp-1:nyp+1)
          yy(1:3)=y(nyp-1:nyp+1)
        endif

        call util_parabel(xx,yy,a,yp,xopt,yopt,ifail)

        if (ifail.eq.0) then

          if (nyp.eq.1) then
            y1(1)=yp(1)
          else if (nyp.eq.n) then
            y1(n)=yp(3)
          else
            y1(nyp)=yp(2)
          endif
          if (n.eq.3) then
            y1(1)=yp(1)
            y1(2)=yp(2)
            y1(3)=yp(3)
            return
          endif
        else
          ifail=-2
          return
        endif
      endif

      do j=nyp,2,-1
        yppl=2.0d0*(y(j-1)-y(j)-y1(j)*(x(j-1)-x(j)))/
     &    (x(j-1)-x(j))**2 ! exchanged
          y1(j-1)=y1(j)+yppl*(x(j-1)-x(j))
      enddo

      do j=nyp,n-1
        yppl=2.0d0*(y(j+1)-y(j)-y1(j)*(x(j+1)-x(j)))/
     &    (x(j+1)-x(j))**2
        y1(j+1)=y1(j)+yppl*(x(j+1)-x(j))
      enddo

      ifail=0

      return
      end
