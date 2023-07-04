*CMZ : 00.00/16 12/09/2014  15.08.59  by  Michael Scheer
*-- Author :    Michael Scheer   12/09/2014
      subroutine util_getarg_int(narg,ival)

      implicit none

      integer narg,ival
      character(64) charg

      call getarg(narg,charg)
      read(charg,*)ival

      return
      end
