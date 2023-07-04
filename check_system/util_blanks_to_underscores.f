*CMZ :          08/03/2018  14.14.52  by  Michael Scheer
*CMZ : 00.00/15 20/12/2013  13.21.33  by  Michael Scheer
*-- Author :    Michael Scheer   20/12/2013
      subroutine util_blanks_to_underscores(chin,chout)

c +PATCH,//UTIL/FOR
c +DECK,util_blanks_to_underscores.

      implicit none

        integer lenin,i,ic

        character(*) chin,chout
        character c1
        equivalence (ic,c1)
        ic=0

        lenin=len_trim(chin)
        chout=chin

        do i=1,lenin
          c1=chin(i:i)
          if (ic.le.32.or.ic.ge.127) chout(i:i)='_'
        enddo

      return
      end
