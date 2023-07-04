*CMZ :          15/05/2019  15.57.36  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_copy_file_py(chsrc,chdest,istat)

      implicit none

      integer luno,istat,nfirst,nlast
      real r(1)

      character(*) chsrc,chdest
      character(1024) chpy,cline
      character(16) chran

      call random_number(r)
      write(chran,*)int(r*10**8)
      call  util_string_trim(chran,nfirst,nlast)

      chpy=".rcopy" // chran(nfirst:nlast) // ".py"
      open(newunit=luno,file=trim(chpy))
      write(luno,'(a)')
     &  'import shutil'
      write(luno,'(a)')
     &  'shutil.copy("'// trim(chsrc) // ',' // trim(chdest) //'");'
      flush(luno)
      close(luno)

      cline="python3 " // trim(chpy)
      istat=system(trim(cline))

      return
      end
