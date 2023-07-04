*CMZ :  4.00/17 28/11/2022  15.15.17  by  Michael Scheer
*CMZ :  4.00/11 21/11/2020  12.01.26  by  Michael Scheer
*CMZ :  4.00/07 12/05/2020  13.45.30  by  Michael Scheer
*CMZ :  4.00/04 19/11/2019  17.28.12  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine wave_check_pid(ipid,istat)

      !use waveenv

      implicit none
      integer ipid,istat,nfirst,nlast,luno,kstat
      real r(1)
      character(16) chran, chpid, cstat
      character(1024) chpy,cline,chpyout

*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.

      write(chpid,*) ipid

      if (kpython.ne.0) then
        call random_number(r)
        write(chran,*)int(r*10**8)
        call  util_string_trim(chran,nfirst,nlast)
        chpy=".rcheckpid" // chran(nfirst:nlast) // ".py"
        chpyout=".rcheckpid" // chran(nfirst:nlast) // ".out"
        write(chpid,*) ipid
        open(newunit=luno,file=trim(chpy))
        write(luno,'(a)') 'import psutil'
        write(luno,'(a)') 'fout = open("' // trim(chpyout) // '","w")'
        write(luno,'(a)') 'fout.write(str(psutil.pid_exists(' // trim(chpid) //'))' // '+ "\n")'
        flush(luno)
        close(luno)
        cline=trim(chpythonhome) // chpathsep // trim(chpythoncom) // " " // trim(chpy)
        istat=system(trim(cline))
        open(newunit=luno,file=trim(chpyout))
        read(luno,*) cstat
        close(luno)
        if (cstat.eq.'True') then
          istat = 1
        else
          istat = 0
        endif
        call util_file_delete(trim(chpy),kstat)
        call util_file_delete(trim(chpyout),kstat)
      else
        if (chplatform.eq.'WINDOWS') then
          call util_check_pid_windows(ipid,istat)
        else if (chplatform.eq.'LINUX' .or. chplatform.eq.'MINGW') then
          call util_check_pid_linux(ipid,istat)
        else
          print*,
     &      "*** Warning in wave_check_pid: No python, don't know how to check process ***"
          istat=-1
        endif
      endif

      return
      end
