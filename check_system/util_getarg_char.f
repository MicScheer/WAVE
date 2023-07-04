*CMZ : 00.00/16 12/09/2014  15.08.59  by  Michael Scheer
*-- Author :    Michael Scheer   12/09/2014
      subroutine util_getarg_real(narg,charg)

      implicit none

      integer narg
      character(*) charg

      call getarg(narg,charg)

      return
      end
