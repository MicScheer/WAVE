*CMZ :          12/05/2017  16.59.45  by  Michael Scheer
*-- Author :    Michael Scheer   12/05/2017
      subroutine util_minimize_fcn(x,f,ifail)

      implicit none

      double precision x,f
      integer ifail

      ifail=0
      f=x**4

      return
      end
