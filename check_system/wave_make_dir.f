*CMZ :  4.00/17 21/11/2022  15.33.16  by  Michael Scheer
*CMZ :  4.00/11 21/11/2020  11.27.00  by  Michael Scheer
*CMZ :  4.00/04 28/06/2019  13.21.55  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine wave_make_dir(chdir,istat)
      !use waveenv

      implicit none

*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.

      integer luno,istat,kstat,nfirst,nlast,iutil_fexist
      real r(1)

      character(*) chdir
      character(1024) chpy,cline
      character(16) chran

      if (iutil_fexist(trim(chdir)).ne.0) then
        istat=1
        return
      endif

      if (kpython.ne.0) then
        call random_number(r)
        write(chran,*)int(r*10**8)
        call  util_string_trim(chran,nfirst,nlast)
        chpy=".mkdir" // chran(nfirst:nlast) // ".py"
        open(newunit=luno,file=trim(chpy))
        write(luno,'(a)')
     &    'import os'
        write(luno,'(a)')
     &    'if os.path.exists("' // trim(chdir) // '") == False:'
        write(luno,'(a)')
     &    '    os.mkdir("' // trim(chdir) // '");'
        flush(luno)
        close(luno)
        cline=trim(chpythonhome) // chpathsep // trim(chpythoncom) // " " // trim(chpy)
        istat=system(trim(cline))
        call util_file_delete(trim(chpy),kstat)
      else
        istat=system(trim(chmkdir) // " " // trim(chdir))
      endif !kpython

      istat=-1
      if (iutil_fexist(trim(chdir)).ne.0) then
        istat=0
      endif

      return
      end
