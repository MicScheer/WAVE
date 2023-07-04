*CMZ :          15/08/2018  14.59.09  by  Michael Scheer
*CMZ : 00.00/16 20/07/2015  09.58.29  by  Michael Scheer
*-- Author :    Michael Scheer   20/07/2015
      subroutine util_solve_2x2(a,x,ifail)

      implicit none

      real*8 a(2,2),ws(2,2),x(2)
      integer ifail

      call deqn(2,a,2,ws,ifail,1,x)

      return
      end
