*CMZ :  4.00/17 28/11/2022  15.02.19  by  Michael Scheer
*CMZ :  4.00/11 21/11/2020  11.38.49  by  Michael Scheer
*CMZ :  4.00/07 27/04/2020  15.32.49  by  Michael Scheer
*CMZ :  4.00/04 28/06/2019  13.30.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine wave_copy_file(chsrc,chdest,istat)
      !use waveenv

      implicit none

*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.

      integer luno,istat,nfirst,nlast,kstat,iutil_fexist
      real r(1)

      character(*) chsrc,chdest
      character(1024) chpy,cline
      character(16) chran

      if (iutil_fexist(chsrc).eq.0) then
        print*,""
        print*,"*** Error in wave_copy_file: File not found!"
        print*,"*** " // trim(chsrc)
        print*,""
        stop "*** Program WAVE aborted ***"
      endif

      if (kpython.ne.0) then

        call random_number(r)
        write(chran,*)int(r*10**8)
        call  util_string_trim(chran,nfirst,nlast)

        chpy=".rcopy" // chran(nfirst:nlast) // ".py"
        open(newunit=luno,file=trim(chpy))
        write(luno,'(a)')
     &    'import shutil'
        write(luno,'(a)')
     &    'shutil.copy("'// trim(chsrc) // '","' // trim(chdest) //'");'
        flush(luno)
        close(luno)

        cline=trim(chpythonhome) // chpathsep // trim(chpythoncom) // " " // trim(chpy)
        istat=system(trim(cline))

        call util_file_delete(trim(chpy),kstat)

      else
        call system(
     &    trim(chcpfile) // " " // trim(chsrc) // " " // trim(chdest))
      endif

      return
      end
