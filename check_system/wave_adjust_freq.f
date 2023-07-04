*CMZ :  3.03/02 11/01/2016  09.42.50  by  Michael Scheer
*-- Author :    Michael Scheer   11/01/2016
*CMZ :          14/12/2015  16.32.39  by  Michael Scheer
*-- Author :    Michael Scheer   14/12/2015
      program wave_adjust

c     Utility to write the current results of wave.adjust to wave_adjust.dat

      implicit none

      double precision chi2,fd0,eopt
      integer i,lenline,nargs

      character(1024) chfile
      character(1024) cline
      character(256) chkey

      chfile="wave.out"

      nargs=command_argument_count()
      if (nargs.eq.1) then
        call get_command_argument(1,chfile)
      endif

      open(unit=20,file=trim(chfile),status='old')
1     read(20,'(a)') cline
      lenline=len_trim(cline) ! length of cline without trailing blanks

      chkey="Estimated maximum:"
      i=index(cline,trim(chkey)) !position of first occurance of chkey in cline
      !read the value from cline after the keyword
      if (i.gt.0) then
        read(20,*) eopt,fd0
c        read(cline(i+len_trim(chkey):lenline),*) fd0
        goto 9
      endif

      goto 1

9     close(20)

c This chi2 will be minimized by WAVE
c      chi2=(1.0d20/fd0)**2
      chi2=(eopt-150000.0d0)**2

      open(unit=20,file='wave_adjust.dat')
      write(20,*)chi2
      close(20)

      open(unit=20,file='wave_adjust.lis',access='append')
      write(20,*)eopt,fd0,chi2

      end
