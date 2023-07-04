*CMZ : 00.00/07 08/05/2008  15.03.41  by  Michael Scheer
*-- Author :    Michael Scheer   08/05/2008
      subroutine util_interpol_parabel_coef(n,xa,ya,a3n,ifail)

      implicit none

      double precision xa(n),ya(n),a3n(3,n),a(3),x(3),y(3),xopt,yopt,yp(3)
      integer i,n,ifail

      do i=1,n-2
        x(1:3)=xa(i:i+2)
        y(1:3)=ya(i:i+2)
        call util_parabel(x,y,a,yp,xopt,yopt,ifail)
        if (ifail.ne.0) then
          return
        endif
        a3n(1:3,i+1)=a(1:3)
      enddo

      a3n(1:3,1)=a3n(1:3,2)
      a3n(1:3,n)=a3n(1:3,n-1)

      return
      end
