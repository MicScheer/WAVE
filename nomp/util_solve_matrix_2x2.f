*CMZ :  3.05/11 15/08/2018  15.03.26  by  Michael Scheer
*-- Author :    Michael Scheer   15/08/2018
*CMZ :          11/01/2018  11.54.22  by  Michael Scheer
*-- Author :    Michael Scheer   10/01/2018
      subroutine util_solve_matrix_2x2(a,b,x,ifail)

c +PATCH,//UTIL/UTIL
c +DECK,util_solve_matrix_2x2.

      ! Calculate matrix x, such that b=x*a

      implicit none

      double precision x(2,2),a(2,2),b(2,2),det
      integer :: ifail

      ifail=0
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if (abs(det).lt.1.0d-30) then
        ifail=-1
        return
      endif

      x(1,1)=(-a(2,1)*b(1,2)+a(2,2)*b(1,1))/det
      x(1,2)=(+a(1,1)*b(1,2)-a(1,2)*b(1,1))/det
      x(2,1)=(-a(2,1)*b(2,2)+a(2,2)*b(2,1))/det
      x(2,2)=(+a(1,1)*b(2,2)-a(1,2)*b(2,1))/det

c      print*,x(1,1)*a(1,1)+x(1,2)*a(2,1)
c      print*,x(1,1)*a(1,2)+x(1,2)*a(2,2)
c      print*,x(2,1)*a(1,1)+x(2,2)*a(2,1)
c      print*,x(2,1)*a(1,2)+x(2,2)*a(2,2)

      return
      end
