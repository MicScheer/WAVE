*CMZ : 00.00/07 19/03/2010  16.41.49  by  Michael Scheer
*CMZ : 00.00/05 07/03/2007  10.40.20  by  Michael Scheer
*-- Author :    Michael Scheer   07/03/2007
      subroutine util_string_firstcharacter(cstring,ifirst,cfirst,istat)

c uses string-handling package M432 from CERN

      implicit none

        integer ilen,ifirst,istat

        character(*) cstring
        character cfirst

*KEEP,strings.
      include 'strings.cmn'
*KEND.

        ilen=len(cstring)

        istat=1
        ifirst = icfnbl(cstring,1,ilen)
        if (ifirst.le.0.and.ifirst.gt.ilen) return

        istat=0
        cfirst=cstring(ifirst:ifirst)

      return
      end
