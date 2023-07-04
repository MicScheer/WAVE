*CMZ : 00.00/15 29/08/2012  11.50.22  by  Michael Scheer
*CMZ : 00.00/07 10/06/2008  14.53.45  by  Michael Scheer
*CMZ : 00.00/02 07/01/98  15.23.41  by  Michael Scheer
*-- Author :    Michael Scheer   07/01/98
      subroutine util_hbook_file_ini(lun,file,cname,istat)

      implicit none

        integer ndpawcp,istat,lrec,lun,iquest,lasth

        external function util_igetlastchar
        integer util_igetlastchar

      parameter (ndpawcp=200000)
        real*4 rpaw(ndpawcp)

        character c
        character(*) file
        character(8) cname

      common/pawc/rpaw
      common/quest/iquest(100)

        call hlimit(ndpawcp)

        iquest(10)=2**17
        lrec=1024

c        lasth=util_igetlastchar(1,1024,file,c)
        lasth=len_trim(file)

        call hropen(lun,cname,file(1:lasth),'nq',lrec,istat)

      return
      end
