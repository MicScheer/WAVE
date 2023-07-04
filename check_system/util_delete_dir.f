*CMZ :          15/05/2019  09.22.24  by  Michael Scheer
*-- Author : Michael Scheer

! Deletes directory chdir using python3
! For kempty.ne.0, directory must be empty

      subroutine util_delete_dir(chdir,kempty,istat)

      implicit none

      integer luno,istat,kstat,nfirst,nlast,iutil_fexist,kempty
      real r(1)

      character(*) chdir
      character(1024) chpy,cline
      character(16) chran

      if (iutil_fexist(trim(chdir)).eq.0) then
        istat=0
        return
      endif

      call random_number(r)
      write(chran,*)int(r*10**8)
      call  util_string_trim(chran,nfirst,nlast)

      if (kempty.ne.0) then
        chpy=".rmemptydir" // chran(nfirst:nlast) // ".py"
        open(newunit=luno,file=trim(chpy))
        write(luno,'(a)')
     &    'import os'
        write(luno,'(a)')
     &    'os.rmdir("' // trim(chdir) // '");'
        flush(luno)
        close(luno)
      else
        chpy=".rmtree" // chran(nfirst:nlast) // ".py"
        open(newunit=luno,file=trim(chpy))
        write(luno,'(a)')
     &    'import shutil'
        write(luno,'(a)')
     &    'shutil.rmtree("' // trim(chdir) // '");'
        flush(luno)
        close(luno)
      endif !kempty

      cline="python3 " // trim(chpy)
      istat=system(trim(cline))

      call util_file_delete(trim(chpy),kstat)

      istat=-1
      if (iutil_fexist(trim(chdir)).eq.0) then
        istat=0
      endif

      return
      end
