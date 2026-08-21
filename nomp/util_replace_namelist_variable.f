*CMZ :          21/08/2026  08.28.08  by  Michael Scheer
*-- Author :    Michael Scheer   20/08/2026
      subroutine util_replace_namelist_variable(cline,chvar,chval,cout,istat)

      implicit none

      integer i,ifound,istat

      character(*) cline,chvar,chval,cout
      character(2048) clinevar,clinecom,chwork
      character c1

      istat=0

      ifound=0

      do i=1,len(cline)
        if (cline(i:i).eq.'!') then
          ifound=i
          exit
        endif
      enddo

      if (ifound.eq.1) then
        istat=1
        return
      endif

      if (ifound.gt.1) then
        clinevar=cline(1:ifound-1)
        clinecom=cline(ifound:len_trim(cline))
      endif

      ifound=0
      do i=1,len(clinevar)
        if (clinevar(i:i).eq.'=') then
          ifound=i
          exit
        endif
      enddo

      if (ifound.eq.0) then
        istat=2
        return
      endif

      chwork=adjustl(clinevar(1:ifound-1))

      if (trim(chwork).eq.chvar) then
        cout=clinevar(1:ifound) // adjustl(trim(chval)) // " " // trim(clinecom)
      endif

c      write(6,'(a)') trim(cout)
      return
      end
