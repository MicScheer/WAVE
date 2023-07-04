*CMZ :  4.00/15 13/03/2022  17.05.32  by  Michael Scheer
*CMZ :  3.03/02 01/11/2016  12.45.16  by  Michael Scheer
*CMZ : 00.00/06 07/01/2008  14.32.30  by  Michael Scheer
*CMZ :  1.19/07 22/08/2002  15.44.21  by  Michael Scheer
*-- Author :    Michael Scheer   09/11/2001
      subroutine util_string_trim(cline,nfirst,nlast)

      implicit none

      integer nfirst,i,ic,nlast
      character(*) cline
      character c1

      equivalence (ic,c1)

      nfirst=-1
      nlast=-1

      do i=1,len(cline)
        c1=cline(i:i)
        if (c1.ne.' '.and.ic.ne.9) then !no blanks, no tabs
          nfirst=i
          goto 1
        endif
      enddo

1     if (nfirst.le.0) return

      do i=len(cline),nfirst,-1
        c1=cline(i:i)
        if (c1.ne.' '.and.ic.ne.9) then !no blanks, no tabs
          nlast=i
          return
        endif
      enddo

      return
      end
