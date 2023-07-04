*CMZ :  4.01/02 07/05/2023  12.16.10  by  Michael Scheer
*CMZ :  4.00/16 08/08/2022  13.42.37  by  Michael Scheer
*CMZ :  4.00/09 10/08/2020  10.21.02  by  Michael Scheer
*CMZ :  4.00/07 06/06/2020  13.57.14  by  Michael Scheer
*CMZ :  3.05/01 03/05/2018  14.17.14  by  Michael Scheer
*CMZ :  3.05/00 02/05/2018  10.55.30  by  Michael Scheer
*-- Author :    Michael Scheer   25/04/2018
      subroutine bundumag_ini(xsta,xsto,lungfo)

      use bpolyederf90m
      use undumagf90m

      implicit none

*KEEP,undumagc.
      include 'undumagc.cmn'
*KEND.

      double precision xsta,xsto,dum,xstao,xstoo

      integer lungfo,lundu,ieof,imag,nplan,ncorn,iplan,icorn,imat,ncol,
     &  n1bpebc,n1bperot,n2bpetm,iw,mag,i

      integer :: idebug=0,istat=0

      character(2048) cline

      logical lexist

      if (ibunduini.eq.1.or.kbundumag_c.ne.0) then
        xstoo=xsto
        xstao=xsta
        if (kbundumag_c.eq.1) then
          perlenclc=4.0d0*(umaglx+uairgap)/1000.0d0
        else if (kbundumag_c.eq.3) then
          perlenclc=uperlen_h/1000.0d0
        endif

        kundumap=kbundumap_c
        call bundumap(0.0d0,0.0d0,0.0d0,dum,dum,dum,xsta,xsto)

        if (xstao.eq.9999.0d0) then
          if (kbundumag_c.eq.1.and.nunduperw.gt.nunduper) then
            xsta=xsta-perlenclc*(nunduperw-nunduper)/2.0d0
          else if (kbundumag_c.eq.3) then
            xsta=xsta-perlenclc*(nperiod_h-nperiodw_h)/2.0d0
          endif
        endif
        if (xstoo.eq.9999.0d0) then
          if (kbundumag_c.eq.1.and.nunduperw.gt.nunduper) then
            xsto=xsto+perlenclc*(nunduperw-nunduper)/2.0d0
          else if (kbundumag_c.eq.3) then
            xsto=xsto+perlenclc*(nperiod_h-nperiodw_h)/2.0d0
          endif
        endif
        return
      endif

      kwave=1

      inquire(file='undumag.wav',exist=lexist)

      if (lexist.eqv..false.) then
        write(lungfo,*)""
        write(lungfo,*)"*** Error in bundumag_ini: File undumag.wav not found ***"
        write(lungfo,*)"Did you run WAVE with KBUNDUMAG=2 before?"
        write(lungfo,*)""
        print*,""
        print*,"*** Error in bundumag_ini: File undumag.wav not found ***"
        print*,"Did you run WAVE with KBUNDUMAG=2 before?"
        print*,""
        stop "*** Program WAVE aborted ***"
      endif

      open(newunit=lundu,file='undumag.wav')

      if (idebug.ne.0) then
        read(lundu,'(a)')cline
        print*,trim(cline)
        backspace(lundu)
      endif

      read(lundu,'(a)')cundutit
      cundutit(1:1)=' '
      read(cundutit,*)kundurun

      if (idebug.ne.0) then
        read(lundu,'(a)')cline
        print*,trim(cline)
        backspace(lundu)
      endif
      read(lundu,*)ncwires,nrec,nmag,nplanmax,ncornmax,
     &  n1bpebc,n1bperot,n2bpetm

      if (idebug.ne.0) then
        read(lundu,'(a)')cline
        print*,trim(cline)
        backspace(lundu)
      endif
      read(lundu,*)ixsym,iysym,izsym,kxcenter,xsym,xcenter

      if (idebug.ne.0) then
        read(lundu,'(a)')cline
        print*,trim(cline)
        backspace(lundu)
      endif

      read(lundu,*,iostat=istat) perlen,xmapmin,xmapmax,nunduclc
      perlenclc=perlen/1000.0d0
      if (istat.ne.0) then
        if (kbundumag_c.eq.1 .or. kbundumag_c.eq.2) then
          nunduclc=nunduper
        else if (kbundumag_c.eq.3 .or. kbundumag_c.eq.4) then
          nunduclc=nperiod_h
        else
          nunduclc=-1
        endif
      endif

      if (idebug.ne.0) then
        read(lundu,'(a)')cline
        print*,trim(cline)
        backspace(lundu)
      endif
      read(lundu,*)kurad,ebeam

      if (idebug.ne.0) then
        read(lundu,'(a)')cline
        print*,trim(cline)
        backspace(lundu)
      endif
      read(lundu,*)xelec,yelec,zelec
      if (idebug.ne.0) then
        read(lundu,'(a)')cline
        print*,trim(cline)
        backspace(lundu)
      endif
      read(lundu,*)vxelec,vyelec,vzelec
      if (idebug.ne.0) then
        read(lundu,'(a)')cline
        print*,trim(cline)
        backspace(lundu)
      endif
      read(lundu,*)xf,yf,zf
      if (idebug.ne.0) then
        read(lundu,'(a)')cline
        print*,trim(cline)
        backspace(lundu)
      endif
      read(lundu,*)efx,efy,efz
      if (idebug.ne.0) then
        read(lundu,'(a)')cline
        print*,trim(cline)
        backspace(lundu)
      endif
      read(lundu,*)tiny,window

      allocate(
     &  bpebc(n1bpebc,nmag),
     &  bpemag(n1bperot,ncornmax,nplanmax,nmag),
     &  bperot(n1bperot,ncornmax,nplanmax,nmag),
     &  bpetm(n1bperot,n2bpetm,nplanmax,nmag),
     &  ibpeplan(nmag),ibpecol(nmag),
     &  ibpecorn(nplanmax,nmag))

      do imag=1,nmag

        if (idebug.ne.0) then
          read(lundu,'(a)')cline
          print*,trim(cline)
          backspace(lundu)
        endif
        read(lundu,*)mag

        if (mag.ne.imag) stop "*** Bad file undumag.wav ***"

        if (idebug.ne.0) then
          read(lundu,'(a)')cline
          print*,trim(cline)
          backspace(lundu)
        endif
        read(lundu,*)bpebc(1:n1bpebc,mag)

        if (idebug.ne.0) then
          read(lundu,'(a)')cline
          print*,trim(cline)
          backspace(lundu)
        endif
        read(lundu,*)nplan,ibpecol(mag)

        ibpeplan(mag)=nplan
        do iplan=1,nplan

          if (idebug.ne.0) then
            read(lundu,'(a)')cline
            print*,trim(cline)
            backspace(lundu)
          endif
          read(lundu,*)ncorn

          ibpecorn(iplan,mag)=ncorn
          do icorn=1,ncorn

            if (idebug.ne.0) then
              read(lundu,'(a)')cline
              print*,trim(cline)
              backspace(lundu)
            endif
            read(lundu,*)bpemag(1:n1bperot,icorn,iplan,mag)

            if (idebug.ne.0) then
              read(lundu,'(a)')cline
              print*,trim(cline)
              backspace(lundu)
            endif
            read(lundu,*)bperot(1:n1bperot,icorn,iplan,mag)

          enddo !ncorn
          do i=1,n2bpetm
            if (idebug.ne.0) then
              read(lundu,'(a)')cline
              print*,trim(cline)
              backspace(lundu)
            endif
            read(lundu,*)bpetm(1:n1bperot,i,iplan,mag)
          enddo
        enddo !nplan
      enddo

      allocate(wire(nwitems,ncwires))

      do iw=1,ncwires
        if (idebug.ne.0) then
          read(lundu,'(a)')cline
          print*,trim(cline)
          backspace(lundu)
        endif
        read(lundu,*) wire(:,iw)
      enddo

      close(lundu)

      niron=nmag-nrec

      write(lungfo,*)
      write(lungfo,*)'     Header of undumag.mag:'
      write(lungfo,*)'     ',trim(cundutit)
      write(lungfo,*)
      write(lungfo,*)'     UNUDMAG variables read from file:'
      write(lungfo,*)
      write(lungfo,*)'     PerLen, xMapMin, xMapMax [mm]:',perlen,xmapmin,xmapmax
      write(lungfo,*)'     kUrad, eBeam:',kurad,ebeam
      write(lungfo,*)'     xElec, yElec, zelec [m]:',xelec,yelec,zelec
      write(lungfo,*)'     VxElec, VyElec, VzElec:',xelec,yelec,zelec
      write(lungfo,*)'     Xf, Yf, Zf:',xf,yf,zf
      write(lungfo,*)'     efX, efY, efZ:',efx,efy,efz
      write(lungfo,*)
      write(lungfo,*)'     Number of magnetic and iron voxels, sum:',
     &  nrec,niron,nmag
      write(lungfo,*)'     Number of current filaments:',ncwires
      write(lungfo,*)
      write(lungfo,*)'     uCorrTiny, tiny, window [mm]:   ',ucorrtiny,tiny,window
      write(lungfo,*)
      write(lungfo,*)'     uRandoX, uRandoY, uRandoZ [mm]:',urandox,urandoy,urandoz
      write(lungfo,*)

      if (xsta.eq.9999.0d0) then
        if (kxstart.eq.1) then
          xsta=xelec
        else if (kxstart.eq.-1) then
          xsta=xmapmin/1000.0d0
        endif
        if (kbundumag_c.eq.1 .or. kbundumag_c.eq.2 .and. nunduperw.gt.nunduclc) then
          xsta=xsta-perlenclc*(nunduperw-nunduclc)/2.0d0
        else if (kbundumag_c.eq.3 .or. kbundumag_c.eq.4 .and. nperiodw_h.gt.nunduclc) then
          xsta=xsta-perlenclc*(nperiodw_h-nunduclc)/2.0d0
        endif
      endif

      if (xsto.eq.9999.0d0) then
        if (kxstop.eq.1) then
          xsto=xf
        else if (kxstop.eq.-1) then
          xsto=xmapmax/1000.0d0
        endif
        if (kbundumag_c.eq.1 .or. kbundumag_c.eq.2 .and. nunduperw.gt.nunduclc) then
          xsto=xsto+perlenclc*(nunduperw-nunduclc)/2.0d0
        else if (kbundumag_c.eq.3 .or. kbundumag_c.eq.4 .and. nperiodw_h.gt.nunduclc) then
          xsto=xsto+perlenclc*(nperiodw_h-nunduclc)/2.0d0
        endif
      endif

      ucorrtinymm=ucorrtiny
      urandoxmm=urandox
      urandoymm=urandoy
      urandozmm=urandoz

      ucorrtiny=ucorrtiny/1000.0d0
      urandox=urandox/1000.0d0
      urandoy=urandoy/1000.0d0
      urandoz=urandoz/1000.0d0

      corrtiny=ucorrtiny

      if (uwwindow.lt.0.0d0) uwindow=1.0d30
      uwindow=uwwindow
      window=uwindow

      ibunduini=1

      return
      end
