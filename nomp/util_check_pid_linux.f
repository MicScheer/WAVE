*CMZ :  4.00/07 12/05/2020  13.53.23  by  Michael Scheer
*CMZ :  4.00/04 19/11/2019  17.28.12  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_check_pid_linux(pid,istat)

      implicit none

      integer pid,istat,nfirst,nlast,luno,kstat,lpid,ifound,i
      real r(1)
      character(16) chran, chpid, cstat
      character(1024) chpy,cline,chpyout

      write(chpid,*) pid
      call  util_string_trim(chpid,nfirst,nlast)
      chpid=chpid(nfirst:nlast)
      lpid = nlast-nfirst+1

      call random_number(r)
      write(chran,*)int(r*10**8)
      call  util_string_trim(chran,nfirst,nlast)
      chpy=".w" // chran(nfirst:nlast) // ".pid"
      istat = system("ps -a > " // trim(chpy))
      if (istat.ne.0) then
        istat = -1
        return
      endif
      open(newunit=luno,file=trim(chpy))
      ifound=0
      do while (.true.)
        read(luno,*,end=91) cstat
        call  util_string_trim(cstat,nfirst,nlast)
        cstat=cstat(nfirst:nlast)
        do i=1,16
          if (cstat(i:i).eq.' ') exit
        enddo
        if (cstat(1:i-1) .eq. chpid(1:lpid)) then
          ifound=1
          exit
        endif
      enddo

 91   close(luno)

      istat = system("rm " // trim(chpy))
      if (istat.ne.0) then
        istat = -1
        return
      endif

      istat=ifound

      return
      end
