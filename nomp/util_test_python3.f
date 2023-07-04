*CMZ :  4.00/17 03/11/2022  12.52.10  by  Michael Scheer
*CMZ :  4.00/07 12/05/2020  09.51.13  by  Michael Scheer
*CMZ :  4.00/04 28/06/2019  12.34.52  by  Michael Scheer
*-- Author : Michael Scheer
! Check python3
      subroutine util_test_python3(chpycom,ifound)

      implicit none

      integer luno,ifound,nfirst,nlast,kstat
      real r(1)

      character(*) chpycom
      character(1024) chpy,cline
      character(16) chran

      if (trim(chpycom).eq.'none') then
        ifound=0
        return
      endif

      call random_number(r)
      write(chran,*)int(r*10**8)
      call  util_string_trim(chran,nfirst,nlast)
      chpy=".testpython3" // chran(nfirst:nlast) // ".py"

      open(newunit=luno,file=trim(chpy))

      write(luno,'(a)') 'import sys'
      write(luno,'(a)') 'sys.exit()'
      flush(luno)
      close(luno)

      if (trim(chpycom).eq.'auto') then
        cline='python3 ' // trim(chpy)
        ifound=system(trim(cline))
        if (ifound.eq.0) then
          ifound=1
          chpycom='python3'
        else
          cline='/usr/bin/python3 ' // trim(chpy)
          ifound=system(trim(cline))
          if (ifound.eq.0) then
            ifound=1
            chpycom='/usr/bin/python3'
          else
            ifound=0
          endif
        endif
      else
        cline=trim(chpycom) // " " // trim(chpy)
        ifound=system(trim(cline))
        if (ifound.eq.0) then
          ifound=1
        else
          ifound=0
        endif
      endif

      call util_file_delete(trim(chpy),kstat)

      return
      end
