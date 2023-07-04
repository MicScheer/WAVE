*CMZ : 00.00/15 20/12/2013  13.09.42  by  Michael Scheer
*-- Author :    Michael Scheer   20/12/2013
      subroutine util_blanks_to_backslash_blanks(chin,chout)

c +PATCH,//UTIL/FOR
c +DECK,util_blanks_to_backslash_blanks.

      implicit none

        integer k,i

        character(*) chin,chout

        k=0

        do i=1,len_trim(chin)
          if (chin(i:i).eq.' ') then
            chout(k+1:k+2)='\ '
            k=k+2
          else
            k=k+1
            chout(k:k)=chin(i:i)
          endif
        enddo

      return
      end
