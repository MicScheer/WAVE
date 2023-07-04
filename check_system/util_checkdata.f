*CMZ : 00.00/07 24/08/2010  11.45.54  by  Michael Scheer
*CMZ : 00.00/02 14/08/2006  13.22.55  by  Michael Scheer
*-- Author :    Michael Scheer   23/01/2004
      subroutine util_checkdata(nlast,cline,last,istat)

      implicit none

      integer nlast,last,istat,idig,jdig,i,itab

      character c1,ctab,cblank
      character(nlast) cline

      equivalence (itab,ctab)

      data cblank/' '/
      data itab/9/

      idig=0
      last=0
      istat=0

      do i=1,nlast

        c1=cline(i:i)

        if (
     &      c1.eq.'0'.or.
     &      c1.eq.'1'.or.
     &      c1.eq.'2'.or.
     &      c1.eq.'3'.or.
     &      c1.eq.'4'.or.
     &      c1.eq.'5'.or.
     &      c1.eq.'6'.or.
     &      c1.eq.'7'.or.
     &      c1.eq.'8'.or.
     &      c1.eq.'9'
     &      ) then
          idig=idig+1
          jdig=1
        else
          jdig=0
        endif

        if (c1.eq.'!'.or.c1.eq.'*'.or.c1.eq.'#'.or.c1.eq.'%') then
          if (idig.eq.0 ) then
            last=0
            istat=1 !comment line
          else
            last=i-1
            istat=0
          endif
          return
        endif

        if (
     &      jdig.eq.0.and.
     &      c1.ne.'.'.and.
     &      c1.ne.'+'.and.
     &      c1.ne.'-'.and.
     &      c1.ne.'e'.and.
     &      c1.ne.'d'.and.
     &      c1.ne.'q'.and.
     &      c1.ne.'E'.and.
     &      c1.ne.'D'.and.
     &      c1.ne.'Q'.and.
     &      c1.ne.','.and.
     &      c1.ne.cblank.and.
     &      c1.ne.ctab
     &      ) then
          istat=-1
          return
        else
          if (c1.eq.'.'.or.jdig.eq.1) last=i
        endif

      enddo

      return
      end
