*CMZ : 00.00/16 12/09/2014  15.05.31  by  Michael Scheer
*-- Author :    Michael Scheer   12/09/2014
      subroutine util_getarg_real(narg,rval)

      implicit none

      integer narg
      real rval
      character(64) charg

      call getarg(narg,charg)
      read(charg,*)rval

      return
      end
