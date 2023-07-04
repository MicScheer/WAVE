*CMZ : 00.00/11 09/03/2011  15.33.07  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_spline_minmax(n,x,y,yp,ypp,xopt,yopt,
     &  ws1,ws2,ws3,ws4,ifail)

      implicit none

      integer n,ifail
      double precision
     &  x(n),y(n),yp(n),ypp(n),xopt,yopt,ws1(n),ws2(n),ws3(n),ws4(n)

      call util_spline_coef_deriv(x,y,n,9999.0d0,9999.0d0,yp,ypp,
     &  ws1,ws2,ws3,ws4)

      call util_spline_zero_derivative(n,x,y,yp,ypp,xopt,yopt,ifail)

      return
      end
