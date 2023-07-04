*CMZ :  4.00/17 28/11/2022  15.15.17  by  Michael Scheer
*CMZ :  4.00/04 28/06/2019  13.31.09  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine wave_move_file(chsrc,chdest,istat)
      !use waveenv

      implicit none

*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.

      integer luno,istat,nfirst,nlast,kstat
      real r(1)

      character(*) chsrc,chdest
      character(1024) chpy,cline
      character(16) chran

      if (kpython.ne.0) then

        call random_number(r)
        write(chran,*)int(r*10**8)
        call  util_string_trim(chran,nfirst,nlast)

        chpy=".rcopy" // chran(nfirst:nlast) // ".py"
        open(newunit=luno,file=trim(chpy))
        write(luno,'(a)')
     &    'import shutil'
        write(luno,'(a)')
     &    'shutil.copystat("'// trim(chsrc) // ',' // trim(chdest) //'");'
        flush(luno)
        close(luno)

        cline=trim(chpythonhome) // " " // trim(chpy)
        istat=system(trim(cline))

        call util_file_delete(trim(chpy),kstat)

      else
        call system(
     &    trim(chmvfile) // " " // trim(chsrc) // " " // trim(chdest))
      endif

      return
      end
