*CMZ : 00.00/07 19/03/2010  16.41.49  by  Michael Scheer
*CMZ : 00.00/05 07/03/2007  10.40.20  by  Michael Scheer
*-- Author :    Michael Scheer   07/03/2007
      subroutine util_string_lastcharacter(cstring,last,clast,istat)

c uses string-handling package M432 from CERN

      implicit none

        integer ilen,last,istat

        character(*) cstring
        character clast

*KEEP,strings.
      include 'strings.cmn'
*KEND.

        ilen=len(cstring)

        istat=1
        last = lnblnk(cstring)
        if (last.le.0.and.last.gt.ilen) return

        istat=0
        clast=cstring(last:last)

      return
      end
