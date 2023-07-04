*CMZ : 00.00/16 27/10/2014  16.54.21  by  Michael Scheer
*-- Author :    Michael Scheer   05/09/2014
      subroutine util_random(n,r)

      integer n
      real :: r(n)

      call random_number(r)

      return
      end
