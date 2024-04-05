*CMZ :          22/01/2018  16.49.42  by  Michael Scheer
*CMZ : 00.00/20 19/08/2016  15.08.19  by  Michael Scheer
*CMZ : 00.00/15 04/01/2013  12.22.07  by  Michael Scheer
*CMZ : 00.00/05 27/02/2007  16.32.04  by  Michael Scheer
*CMZ : 00.00/02 04/08/2006  14.56.41  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.58  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_zeit_kommentar(lun,comment)

c to determine date and time and write it to logical unit lun

      implicit none

      integer lun

      character(*) comment

      character spacer(50)
      character(10) dtday,dttime,dtzone
      integer idatetime(8),ilast

      data spacer/50*' '/

      call date_and_time(dtday,dttime,dtzone,idatetime)

      ilast=len_trim(comment)

      write(lun,*)
      write(lun,*)
      if (ilast.gt.0) then
        write(lun,*)comment(1:ilast),spacer,dttime(1:2),':',dttime(3:4),':',dttime(5:6),' '
     &    ,dtday(7:8),'.',dtday(5:6),'.',dtday(3:4)
      else
        write(lun,*)spacer,dttime(1:2),':',dttime(3:4),':',dttime(5:6),' '
     &    ,dtday(7:8),'.',dtday(5:6),'.',dtday(3:4)
      endif
      write(lun,*)

      return
      end
