*CMZ :          08/03/2018  14.14.52  by  Michael Scheer
*CMZ : 00.00/06 10/01/2008  14.23.25  by  Michael Scheer
*-- Author :    Michael Scheer   10/01/2008
      subroutine util_string_contains_letter(cline,istat)

c check, if cline contains letter, i.e. A-Z,a-z

      character(*) cline
      character c1
      integer i1,i,istat
      equivalence(c1,i1)

      i1=0
      istat=-1 !error

      do i=1,len(cline)
        c1=cline(i:i)
        if (i1.ge.65.and.i1.le.90.or.i1.ge.97.and.i1.le.122) then
          istat=1
          return
        endif
      enddo

      istat=0

      return
      end
