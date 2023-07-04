*CMZ : 00.00/15 29/08/2012  11.50.48  by  Michael Scheer
*CMZ : 00.00/07 10/06/2008  14.08.14  by  Michael Scheer
*CMZ : 00.00/02 07/01/98  15.23.41  by  Michael Scheer
*-- Author :    Michael Scheer   07/01/98
      subroutine util_hbook_file_end(lun,cname)

      implicit none

        integer lun,lastn

        external function util_igetlastchar
        integer util_igetlastchar

        character c
        character(*) cname

c        lastn=util_igetlastchar(1,1024,cname,c)
        lastn=len_trim(cname)

        call hrend(cname(1:lastn))
        close(lun)

      return
      end
