*CMZ : 00.00/16 20/07/2015  09.58.29  by  Michael Scheer
*-- Author :    Michael Scheer   20/07/2015
      subroutine util_solve_6x6(a,x,ifail)

      implicit none

      real*8 a(6,6),ws(6,6),x(6)
      integer ifail

      call deqn(6,a,6,ws,ifail,1,x)

      return
      end
