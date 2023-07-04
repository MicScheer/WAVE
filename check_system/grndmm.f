*CMZ :  3.02/03 27/10/2014  11.48.57  by  Michael Scheer
*-- Author :    Michael Scheer   24/10/2014
      subroutine grndmm(rvec,ndim)

      implicit none

      integer ndim
      real rvec(*)

      call util_random(ndim,rvec)

      return
      end
