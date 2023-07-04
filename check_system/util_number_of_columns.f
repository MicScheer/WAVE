*CMZ :          14/03/2017  21.02.42  by  Michael Scheer
*CMZ : 00.00/19 13/08/2015  13.33.09  by  Michael Scheer
*CMZ : 00.00/07 15/10/2010  09.11.16  by  Michael Scheer
*-- Author :    Michael Scheer   15/10/2010
c +PATCH,//UTIL/FOR
c +DECK,util_getncol.
      subroutine util_number_of_columns(file,ncol,istat)

      implicit none

      integer ncol,lun,last
      integer :: ndim=2048, nwords=0, ipos=0, istat
      character(2048) cline
      character(*) file
      double precision x(2048)

      ncol=0
      istat=-1
      open(newunit=lun,file=file,status='old')
      call  util_read_line(lun,cline,last)
      call util_string_split(cline,ndim,nwords,ipos,istat)
      if (last.ne.0) then
        read(cline,*,err=9999)x(1:nwords)
        ncol=nwords
        istat=0
        goto 9999
      endif

9999  close(lun)
      return
      end
