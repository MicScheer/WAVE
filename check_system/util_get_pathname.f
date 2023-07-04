*CMZ : 00.00/16 10/04/2014  09.47.44  by  Michael Scheer
*-- Author :    Michael Scheer   10/04/2014
      subroutine util_get_pathname(cline,path)

c +PATCH,//UTIL/FOR
c +DECK,util_get_pathname.

      implicit none

      integer i
      character(*) cline, path

      do i=len_trim(cline),1,-1
        if (cline(i:i).eq.'/') goto 1
      enddo

1     continue

      path=cline(1:i)

      return
      end
