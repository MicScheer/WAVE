*CMZ :  4.00/06 05/12/2019  13.07.35  by  Michael Scheer
*CMZ : 00.00/16 19/03/2014  12.14.29  by  Michael Scheer
*CMZ : 00.00/15 03/09/2012  09.29.49  by  Michael Scheer
*CMZ : 00.00/06 08/03/2007  14.02.27  by  Michael Scheer
*CMZ : 00.00/05 07/03/2007  12.58.44  by  Michael Scheer
*-- Author :    Michael Scheer   07/03/2007
      subroutine util_string_substring_igncase(cline,substring,ianf,iend,istat)

c Input:
c      cline, substring

c If substring is passed as variable, full length of substring is checked,
c i.e. pending invisible characters are tested as well

c Output:
c      ianf,iend: start and end position of substring, 0 if not found
c      istat: error, i.e. string not found

c Evtl. besser FORTRAN-functions scan oder index benutzen

      implicit none

      integer ilenl,ilens,istat,ianf,iend,i

      character(*) cline,substring
      character (len = :), allocatable ::  clow,sublow

      istat=-1
      ianf=0
      iend=0

      ilenl=len(cline)
      ilens=len(substring)

      clow=cline
      sublow=substring
      call util_lower_case(clow)
      call util_lower_case(sublow)

      if (ilens.gt.ilenl) return

      do i=1,ilenl-ilens+1
        if (clow(i:i+ilens-1).eq.sublow) then
          ianf=i
          iend=ianf+ilens-1
          istat=0
          return
        endif
      enddo

      return
      end
