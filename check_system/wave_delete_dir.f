*CMZ :  4.00/17 21/11/2022  15.32.48  by  Michael Scheer
*CMZ :  4.00/11 21/11/2020  11.56.57  by  Michael Scheer
*CMZ :  4.00/04 28/06/2019  13.21.55  by  Michael Scheer
*-- Author : Michael Scheer

! For kempty.ne.0, directory must be empty

      subroutine wave_delete_dir(chdir,kempty,istat)
      !use waveenv

      implicit none

*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.

      integer luno,istat,kstat,nfirst,nlast,iutil_fexist,kempty
      real r(1)

      character(*) chdir
      character(1024) chpy,cline
      character(16) chran

      if (iutil_fexist(trim(chdir)).eq.0) then
        istat=0
        return
      endif

      if (kpython.ne.0) then
        call random_number(r)
        write(chran,*)int(r*10**8)
        call  util_string_trim(chran,nfirst,nlast)
        if (kempty.ne.0) then
          chpy=".rmemptydir" // chran(nfirst:nlast) // ".py"
          open(newunit=luno,file=trim(chpy))
          write(luno,'(a)')
     &      'import os'
          write(luno,'(a)')
     &      'os.rmdir("' // trim(chdir) // '");'
          flush(luno)
          close(luno)
        else
          chpy=".rmtree" // chran(nfirst:nlast) // ".py"
          open(newunit=luno,file=trim(chpy))
          write(luno,'(a)')
     &      'import shutil'
          write(luno,'(a)')
     &      'shutil.rmtree("' // trim(chdir) // '");'
          flush(luno)
          close(luno)
        endif !kempty

        cline=trim(chpythonhome) // chpathsep // trim(chpythoncom) // " " // trim(chpy)
        istat=system(trim(cline))

        call util_file_delete(trim(chpy),kstat)

      else
        if (kempty.eq.0) then
          istat=system(trim(chrmtree) // " " // trim(chdir))
        else
          istat=system(trim(chrmdir) // " " // trim(chdir))
        endif
      endif

        istat=-1
        if (iutil_fexist(trim(chdir)).eq.0) then
        istat=0
      endif

      return
      end
