*CMZ :  3.03/02 16/02/2017  12.20.02  by  Michael Scheer
*-- Author :    Michael Scheer   16/02/2017
      subroutine util_parabola_to_zero(x,y,a,x0,ifail)

      ! Calculates parabola y(i)=a*(x(i)-x0)**2

      implicit none

      double precision x(2),y(2),a,x0,sq
      integer ifail

      ifail=-1

      if (y(1)*y(2).le.0.0d0.or.y(1).eq.y(2).or.x(1).eq.x(2)) then
        return
      endif

      sq=sqrt(y(1)/y(2))
      x0=(x(1)-x(2)*sq)/(1.0d0-sq)
      a=y(1)/(x(1)-x0)**2

      ifail=0

      return
      end
