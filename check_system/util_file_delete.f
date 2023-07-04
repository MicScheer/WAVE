*CMZ :          01/08/2018  09.18.53  by  Michael Scheer
*CMZ : 00.00/16 23/06/2014  15.51.32  by  Michael Scheer
*CMZ : 00.00/07 21/07/2009  14.58.29  by  Michael Scheer
*CMZ : 00.00/06 12/07/2007  15.45.32  by  Michael Scheer
*-- Author :    Michael Scheer   12/07/2007
      subroutine util_file_delete(file,istat)

      integer istat,lun

      character(*) file
      logical lexist

      istat=-1

      inquire(file=file,exist=lexist)

      if (lexist.eqv..false.) then
        istat=1
        return
      endif

      open(newunit=lun,file=file,status='old',iostat=istat)
      if (istat.eq.0) close(lun,status='delete')

      return
      end
