*CMZ :  4.01/05 16/04/2024  12.38.45  by  Michael Scheer
*CMZ :  2.05/04 16/12/2023  12.06.33  by  Michael Scheer
*CMZ :  2.04/19 16/09/2023  16.33.50  by  Michael Scheer
*CMZ :  2.03/00 26/07/2022  07.55.50  by  Michael Scheer
*CMZ :  2.02/02 01/07/2022  17.30.28  by  Michael Scheer
*CMZ :  2.02/00 29/03/2021  09.26.44  by  Michael Scheer
*CMZ :  1.00/00 19/08/2016  18.24.11  by  Michael Scheer
*CMZ : 00.00/15 04/01/2013  12.22.07  by  Michael Scheer
*CMZ : 00.00/05 27/02/2007  16.32.04  by  Michael Scheer
*CMZ : 00.00/02 04/08/2006  14.56.41  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.58  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_zeit_kommentar_delta(lun,comment,iset)

c to determine date and time and write it to logical unit lun

      implicit none

      integer lun,ilast

      character(*) comment

      character spacer(30)
      character(10) dtday,dttime,dtzone
      character(32) c32
      character(2048) :: cblank=''
      character(2048) cline

      real :: seconds=0.0,secondso=0.0
      integer idatetime(8),iyear,imonth,iday,ihour,iminute,isec
      integer iyearo,imontho,idayo,ihouro,iminuteo,iseco,kday,khour,kminit,ksec

      integer :: ical=0,linlen,iset,iseconds,isecondso

      data spacer/30*' '/

      save

      call date_and_time(dtday,dttime,dtzone,idatetime)

      iyear=idatetime(1)
      imonth=idatetime(2)
      iday=idatetime(3)
      ihour=idatetime(5)
      iminute=idatetime(6)
      isec=idatetime(7)

      ilast=len_trim(comment)

      write(lun,*)
      if (ilast.gt.0) then
        write(cline,*)comment(1:ilast),spacer,dttime(1:2),':',dttime(3:4),':',
     &    dttime(5:6),' ',dtday(7:8),'.',dtday(5:6),'.',dtday(3:4)
      else
        write(cline,*)spacer,dttime(1:2),':',dttime(3:4),':',dttime(5:6),' '
     &    ,dtday(7:8),'.',dtday(5:6),'.',dtday(3:4)
      endif

      ilast=len_trim(cline)

      if (ilast.lt.64) then
        write(lun,'(a)') cblank(1:64-ilast) // cline(1:ilast)
      else
        write(lun,'(a)') trim(cline)
      endif

      if (ical.gt.0.and.iset.eq.0) then
        cline='  delta time: '
        seconds=secnds(secondso)
        iseconds=nint(seconds)
        kday=0
        khour=0
        kminit=0
        ksec=0
        if (iseconds.ge.3600*24) then
          kday=iseconds/(3600*24)
          iseconds=iseconds-kday*3600*24
          write(c32,*)kday
          cline=trim(cline) // trim(c32) // ' days '
        endif
        if (iseconds.ge.3600) then
          khour=iseconds/3600
          iseconds=iseconds-khour*3600
          write(c32,*)khour
          cline=trim(cline) // trim(c32) // ' hours '
        endif
        if (iseconds.ge.60) then
          kminit=iseconds/60
          iseconds=iseconds-kminit*60
          write(c32,*)kminit
          cline=trim(cline) // trim(c32) // ' minutes '
        endif
        write(c32,*)iseconds
        cline=trim(cline) // trim(c32) // ' seconds'
        write(lun,*) trim(cline)
      endif

      write(lun,*)''

      if (ical.eq.0.or.iset.ne.0) then
        secondso=secnds(0.0)
        isecondso=nint(secondso)
        iyearo=iyear
        imontho=imonth
        idayo=iday
        ihouro=ihour
        iminuteo=iminute
        iseco=isec
      endif

      ical=1

      return
      end
