*CMZ :  4.01/04 28/12/2023  13.26.19  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  15.12.26  by  Michael Scheer
*CMZ :  4.01/00 11/02/2023  16.38.29  by  Michael Scheer
*CMZ :  4.00/17 05/12/2022  10.30.41  by  Michael Scheer
*CMZ :  4.00/16 17/09/2022  15.46.32  by  Michael Scheer
*CMZ :  4.00/15 02/06/2022  09.45.10  by  Michael Scheer
*CMZ :  4.00/11 28/06/2021  10.33.06  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine amprep_omp

! Interface WAVE -> urad_phase

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.

      use sourcef90
      use observf90
      use afreqf90
      use bunchmod
      use uradphasemod
      use omp_lib

      implicit none

*KEEP,datetime.
      include 'datetime.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,track0.
      include 'track0.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,observ.
      include 'observ.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,ellip.
      include 'ellip.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,optic.
      include 'optic.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEND.

      double precision step,sigzp,sigz,sigyp,sigy,perlen,shift,ephmin,ephmax,
     &  ebeam,disp,dispp,curr,beffh,beffv,pherror,pw,ph,xbeta,
     &  beth,alphh,betv,alphv,de,disph,dispph,dispv,disppv

      integer :: modewave=1,icohere,nelec,noranone,modebunch,npinz,npiny,
     &  modepin,modesphere,nper,nepho

c      if (user(2).ne.0) then
c        modewave=user(2)
c        print*,"modewave = USER(2)",modewave
c      else
        modewave=1
c      endif

      step=1.0d0/dble(myinum)*1000.0d0

      perlen=phrperl*1000.0d0
      shift=phrshift*1000.0d0
      beffv=phrb0v
      beffh=phrb0h
      nper=kampli

      nepho=nfreq
      ephmin=freq(1)
      ephmax=freq(nfreq)

      if (rpinsph.gt.0.0d0) modesphere=1

      if (ipin.ne.0) then
        npiny=nobsvy
        npinz=nobsvz
      else
        npiny=1
        npinz=1
      endif

      if (ipin.eq.3) then
        ph=pinh*1000.0d0
        pw=pinw*1000.0d0
      else
        ph=(obsv(2,nobsv)-obsv(2,1))*1000.0d0
        pw=(obsv(3,nobsv)-obsv(3,1))*1000.0d0
      endif

      xbeta=phrxbeta
      beth=phrbetah
      alphh=phralphah
      disph=phrdisph
      dispph=phrdispph

      betv=phrbetav
      alphv=phralphav
      dispv=phrdispv
      disppv=phrdisppv

      de=phrespread
      pherror=phrerror

      bunchlen=phrbunlen
      modebunch=iubunch

      if (ibunch.eq.0) then

        icohere=0
        nelec=1
        noranone=1
        modepin=0

        if (ifold.ne.0) then
          if (ifold.ne.1) then
            print*,''
            print*,'--- Warning in ampre_omp: IFOLD set, but not 1, setting it to 1 ---'
            write(lungfo,*)''
            write(lungfo,*)'--- Warning in ampre_omp: IFOLD set, but not 1, setting it to 1 ---'
            ifold=1
          endif

        endif

      else

        if (nbunch.gt.1.and.neinbunch.gt.1) then
          print*,'*** Error in ampre_omp: Both, NBUNCH and NEINBUNCH > 1. This is not allowed here.'
          print*
          write(lungfo,*)'*** Error in ampre_omp: Both, NBUNCH and NEINBUNCH > 1. This is not allowed here.'
          stop "*** Program WAVE terminated ***"
        endif

        if (neinbunch.gt.1) then
          icohere=1
        else
          icohere=0
        endif

        nelec=max(nbunch,nelec)

        if (ibunch.eq.-1) then
          noranone=1
        else
          noranone=0
        endif

        if (ipin.eq.0.or.ipin.eq.1) then
          modepin=0
        else
          modepin=1
        endif

      endif !ibunch

      call urad_phase(
     &  mthreads,nelec,noranone,icohere,modebunch,bunchlen,bunchcharge,ihbunch,
     &  perlen,shift,nper,beffv,beffh,
     &  dmyenergy,dmycur,step,nlpoi,
     &  pincen*1000.0d0,pw,ph,npiny,npinz,modepin,modesphere,
     &  nepho,ephmin,ephmax,banwid,
     &  xbeta,beth,alphh,betv,alphv,de,phremith,phremitv,
     &  disph,dispph,dispv,disppv,
     &  modeph,pherror,modewave
     &  )

      return
      end
