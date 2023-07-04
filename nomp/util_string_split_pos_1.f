*CMZ :  4.00/07 04/06/2020  09.56.14  by  Michael Scheer
*CMZ : 00.00/06 07/03/2007  17.00.51  by  Michael Scheer
*CMZ : 00.00/05 07/03/2007  12.58.44  by  Michael Scheer
*-- Author :    Michael Scheer   07/03/2007
      subroutine util_string_split_pos_1(cline,ipos,chsep,istat)

c Input:
c      cline, chsep

c Output:
c      ipos: First position of chsep in cline

      implicit none

      integer istat,ipos,jx

      character(*) cline
      character chsep

      istat=-1
      ipos=-1

      do jx=1,len_trim(cline)
        if (cline(jx:jx).eq.chsep) then
          ipos=jx
          istat=0
          exit
        endif
      enddo

      return
      end
