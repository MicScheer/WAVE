*CMZ :  3.05/15 06/10/2018  11.25.38  by  Michael Scheer
*CMZ :  3.05/04 28/06/2018  15.52.47  by  Michael Scheer
*CMZ : 00.00/16 19/03/2014  12.14.18  by  Michael Scheer
*CMZ : 00.00/15 03/09/2012  09.26.58  by  Michael Scheer
*CMZ : 00.00/07 05/03/2008  15.43.44  by  Michael Scheer
*CMZ : 00.00/02 14/08/2006  13.22.55  by  Michael Scheer
*-- Author :    Michael Scheer   23/01/2004
      subroutine util_skip_comment_end(lun,ieof)

      implicit none

      integer lun,ieof
c      character com
      character(32) c2

      ieof=0

1     read(lun,'(a)',end=99) c2

      if (
     &    c2(1:1).ne.'!'.and.c2(1:1).ne.'*'.and.c2(1:1).ne.'#'
     &    .and.c2(1:1).ne.'%'.and.c2(1:1).ne.'@'.and.
     &    c2.ne.' !'.and.c2.ne.' *'.and.c2.ne.' #'
     &    .and.c2.ne.' %'.and.c2.ne.' @'.and.len_trim(c2).ne.0
     &    ) then
        backspace(lun)
      else
        goto 1
      endif

      return

99    ieof=1

      return
      end
