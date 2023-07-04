*CMZ : 00.00/15 19/04/2013  14.45.07  by  Michael Scheer
*-- Author :    Michael Scheer   19/12/2012
      subroutine util_get_hostname(lun,chost,istat)

      integer lun,istat
      character(*) chost

      istat=-1
      chost='unkown host'
      open(unit=lun,file='/etc/hostname',err=9999,status='old')
      read(lun,'(a)')chost
      close(lun)

9999  return

      end
