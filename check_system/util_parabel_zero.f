*CMZ : 00.00/14 17/09/2011  19.20.50  by  Michael Scheer
*CMZ : 00.00/07 22/03/2010  15.28.00  by  Michael Scheer
*CMZ : 00.00/02 26/03/97  10.23.11  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.40  by  Michael Scheer
*-- Author :
      subroutine util_parabel_zero(x,y,a,yp,xopt,yopt,xzero,ifail)

c calculates x0 for 0=a0+a1*(x-x0)+a2*(x-x0)**2

      implicit none

      integer ifail
      double precision x(3),y(3),a(3),yp(3),xzero(2),xopt,yopt,root

      ifail=0

      call util_parabel(x,y,a,yp,xopt,yopt,ifail)

      if (ifail.ne.0) return

      root=-4.0d0*a(1)*a(3)+a(2)**2
      if (root.ge.0) then
        root=sqrt(root)
      else
        ifail=2
      endif
      xzero(1)=(root-a(2))/(2.0d0*a(3))
      xzero(2)=(-root-a(2))/(2.0d0*a(3))

      return
      end
