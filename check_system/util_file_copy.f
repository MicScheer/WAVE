*CMZ : 00.00/07 21/07/2009  14.59.25  by  Michael Scheer
*CMZ : 00.00/06 12/07/2007  15.43.22  by  Michael Scheer
*-- Author :    Michael Scheer   12/07/2007
      subroutine util_file_copy(file1,file2,istat)

      integer istat,len1,len2

      character(*) file1,file2
      character(1000) command

      logical lexist

      istat=-1

      len1=len(file1)
      len2=len(file2)

      if (len1.gt.900.or.len2.gt.900) then
        istat=1
        return
      endif

      inquire(file=file1,exist=lexist)

      if (lexist.eqv..false.) then
        istat=1
        return
      endif


      command='cp '//file1(1:len1)//' '//file2(1:len2)
      print*,command
      call system(command)


      istat=0

      return
      end
