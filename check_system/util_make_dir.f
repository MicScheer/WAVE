*CMZ :          15/05/2019  09.24.19  by  Michael Scheer
*-- Author : Michael Scheer
! Creates directory chdir using python3

      subroutine util_make_dir(chdir,istat)

      implicit none

      integer luno,istat,kstat,nfirst,nlast,iutil_fexist
      real r(1)

      character(*) chdir
      character(1024) chpy,cline
      character(16) chran

      if (iutil_fexist(trim(chdir)).ne.0) then
        istat=1
        return
      endif

      call random_number(r)
      write(chran,*)int(r*10**8)
      call  util_string_trim(chran,nfirst,nlast)
      chpy=".mkdir" // chran(nfirst:nlast) // ".py"

      open(newunit=luno,file=trim(chpy))

      write(luno,'(a)')
     &  'import os'
      write(luno,'(a)')
     &  'if os.path.exists("' // trim(chdir) // '") == False:'
      write(luno,'(a)')
     &  '    os.mkdir("' // trim(chdir) // '");'
      flush(luno)
      close(luno)

      cline="python3 " // trim(chpy)
      istat=system(trim(cline))

      call util_file_delete(trim(chpy),kstat)

      istat=-1
      if (iutil_fexist(trim(chdir)).ne.0) then
        istat=0
      endif

      return
      end
