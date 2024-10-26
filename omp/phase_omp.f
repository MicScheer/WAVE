*CMZ :          20/08/2024  17.20.56  by  Michael Scheer
*CMZ :  4.01/05 11/03/2024  18.40.00  by  Michael Scheer
*CMZ :  4.01/04 28/11/2023  14.20.34  by  Michael Scheer
*CMZ :  4.01/03 12/06/2023  10.59.52  by  Michael Scheer
*CMZ :  4.00/14 07/02/2022  16.17.00  by  Michael Scheer
*CMZ :  3.02/05 22/03/2015  19.55.19  by  Michael Scheer
*CMZ :  3.02/03 23/10/2014  13.43.13  by  Michael Scheer
*CMZ :  3.02/00 15/10/2014  09.29.12  by  Michael Scheer
*CMZ :  3.01/05 12/06/2014  08.52.10  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.13.36  by  Michael Scheer
*CMZ :  2.70/09 15/01/2013  14.58.32  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  15.34.39  by  Michael Scheer
*CMZ :  2.70/03 14/12/2012  14.29.48  by  Michael Scheer
*CMZ :  2.70/02 14/12/2012  10.34.16  by  Michael Scheer
*CMZ :  2.70/01 12/12/2012  15.50.01  by  Michael Scheer
*CMZ :  2.70/00 11/12/2012  17.05.31  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  12.06.21  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  09.55.57  by  Michael Scheer
*CMZ :  2.63/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.49/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.48/04 12/03/2004  15.40.31  by  Michael Scheer
*CMZ :  2.47/21 03/12/2003  09.40.09  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.41/07 13/06/2002  15.07.39  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.20/01 03/01/2001  13.40.19  by  Michael Scheer
*CMZ :  2.16/08 27/10/2000  14.30.15  by  Michael Scheer
*CMZ :  2.16/04 24/06/2000  17.20.05  by  Michael Scheer
*CMZ :  2.16/00 07/06/2000  23.23.42  by  Michael Scheer
*CMZ :  2.15/00 02/05/2000  18.08.47  by  Michael Scheer
*CMZ :  2.13/07 10/02/2000  16.43.36  by  Michael Scheer
*CMZ :  2.13/04 21/01/2000  12.37.13  by  Michael Scheer
*CMZ :  2.13/03 18/01/2000  18.06.22  by  Michael Scheer
*CMZ :  2.13/00 01/12/99  17.14.49  by  Michael Scheer
*CMZ :  2.10/01 30/04/99  13.59.54  by  Michael Scheer
*CMZ :  1.04/00 27/11/98  12.44.48  by  Michael Scheer
*CMZ :  1.03/06 29/09/98  14.43.55  by  Michael Scheer
*-- Author :    Michael Scheer   18/09/98
      SUBROUTINE PHASE_omp
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,phasef90u.
      include 'phasef90u.cmn'
*KEEP,phasewsf90u.
      include 'phasewsf90u.cmn'
*KEEP,wbetaf90u.
      include 'wbetaf90u.cmn'
*KEEP,wbetaf90u.
      include 'wbetaf90u.cmn'
*KEND.

      use ompmod
      use omp_lib
      use wobsvmod

* ROUTINE TO PROPAGATE COMPLEXE AMPLITUDE FROM PINHOLE BACK TO (PHCENX,PHCENY,PHCENZ)


      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,reargf90.
      include 'reargf90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,phasef90.
      include 'phasef90.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,wbetaf90.
      include 'wbetaf90.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,uservar.
      include 'uservar.cmn'
*KEND.

      CHARACTER(8) OLDDIR

      integer :: idebug=1

      INTEGER ICYCLE,NTUP_P,IOBS,IPHZ,iy,iz,iy1,iz1,IPHY,ifrq,IEPS,I,NGEO_P,ISOUR
      INTEGER NIDGEO1,ISTAT,NIDGEO2,NBEAM_P,J,IELEM,NSIZE_P
      INTEGER IOBSY,IOBSZ,K,ix,is,mthreadso

      complex*16 efc(3),bfc(3),expsh,rea(3)

      integer
     &  mphasey_omp,mphasez_omp,nphasey_omp,nphasez_omp,iphfold_omp,
     &  nphelem_omp,ihsel_omp,ith,nfreq_omp,iphase_omp
      integer ic,kobs

      double precision reanor,dum,
     &  phaperzm_omp,phaperzp_omp,phaperzpm_omp,phaperzpp_omp,
     &  phaperym_omp,phaperyp_omp,phaperypm_omp,phaperypp_omp,
     &  phelem_omp(5,4,nphelemp),dgsigz_omp,dgsigy_omp,r(3),ef(3) !,bf(3)

      PARAMETER(NTUP_P=24,NGEO_P=16,NBEAM_P=16,NSIZE_P=4)

      CHARACTER(5) CHTAGS(NTUP_P),CHGEO(NGEO_P),CHBEAM(NBEAM_P)
      CHARACTER(5) CHSIZE(NSIZE_P)

      REAL*8 TUP(NTUP_P),TGEO(NGEO_P),SELGEO,TBEAM(NGEO_P),TSIZ(NSIZE_P)
      REAL*4 FLOW,FHIG,DF

      DOUBLE PRECISION XPH,YPH,ZPH,XOBS,YOBS,ZOBS,DX,DY,DZ,DZY2,ANS
     &  ,OMC,DOMC,DR,DRRED,DX2,DMASHZ,DMASHY,PHLOWZ,PHLOWY,EPS(6)
     &  ,FOCUS,RLAMBDA1,smax,sfmax

      DOUBLE PRECISION XSOUR,YSOUR,ZSOUR,DR2PH,DR2SOUR,THETA,PHI,TANTHE,TANPHI
     &  ,DXPH,DA,EPSBEAM
     &  ,XBEAM,YBEAM,ZBEAM,TANTHEB,TANPHIB,OPTMAT(4,4),BEAM(4)
     &  ,TOTMAT(4,4),DUMMAT(4,4)

      DOUBLE PRECISION W
      DOUBLE PRECISION XA(NDOBSVZP)
     &  ,YAR1(NDOBSVZP)
     &  ,YAR2(NDOBSVZP)
     &  ,YAR3(NDOBSVZP)
     &  ,YAI1(NDOBSVZP)
     &  ,YAI2(NDOBSVZP)
     &  ,YAI3(NDOBSVZP)
     &  ,RESULT(6,2)

      DOUBLE PRECISION xay(ndobsvzp)
     &  ,YAR4(NDOBSVZP)
     &  ,YAR5(NDOBSVZP)
     &  ,YAR6(NDOBSVZP)
     &  ,YAI4(NDOBSVZP)
     &  ,YAI5(NDOBSVZP)
     &  ,YAI6(NDOBSVZP)
     &  ,RESULTY(6,2)

      DOUBLE PRECISION
     &  YAR1Y(NDOBSVYP)
     &  ,YAR2Y(NDOBSVYP)
     &  ,YAR3Y(NDOBSVYP)
     &  ,YAI1Y(NDOBSVYP)
     &  ,YAI2Y(NDOBSVYP)
     &  ,YAI3Y(NDOBSVYP)
     &  ,YAR4Y(NDOBSVYP)
     &  ,YAR5Y(NDOBSVYP)
     &  ,YAR6Y(NDOBSVYP)
     &  ,YAI4Y(NDOBSVYP)
     &  ,YAI5Y(NDOBSVYP)
     &  ,YAI6Y(NDOBSVYP)

      double precision, dimension (:,:,:), allocatable :: phspec3,phspec3f,
     &  phspec3fy

      double precision, dimension (:), allocatable :: zphw,yphw,freq_omp,
     &  specwz,specwy,specfwz,specfwy,phws1,phws2,phcoef,phws3,phws4

      double precision phgsigz,phgsigy,wlen,sigrp,sigr,rn,rnx,rny,rnz

      integer is0

      character(8) chphase
      integer lenchphase

      DATA CHTAGS
     &  /'x','y','z','e','ie','iy','iz',
     &  're_x','im_x','re_y','im_y','re_z','im_z','spec','specf',
     &  'rb_x','ib_x','rb_y','ib_y','rb_z','ib_z',
     &  'nx','ny','nz'
     &  /

      data chgeo
     &  /'x','y','z','yp','zp','e','ie','is','xs','ys','zs','spec'
     &  ,'xo','yo','zo','speco'/
      data chbeam
     &  /'x','y','z','yp','zp','e','ie','is','xs','ys','zs','spec'
     &  ,'xo','yo','zo','speco'/
      data chsize /'ie','is','zrms','yrms'/

      DATA EPSBEAM/0.001D0/

      call zeit(6)
      write(6,*)"     Performing field propagation"

c      print*,"************** kein OMP ************"
c      print*,"*** Vorzeichen von Imag(B) noch korrekt??"

      mthreadso=mthreads
      mthreads=max(1,mthreads)

      if (mhbookp.eq.0) then
        chphase='//PHASE'
        lenchphase=7
      else
        chphase='//WAVE'
        lenchphase=6
      endif

      if (iemit.eq.1) then
        do ix=1,nco-1
          if (wbeta(1,ix).le.phcenx.and.wbeta(1,ix+1).gt.phcenx) then
            is0=ix
          endif
        enddo
      endif

      if (iphfold.ne.0) then
        if (phbeth.eq.-9999.0d0) then
          phbeth=wbeta(2,is0)
        endif
        if (phbetv.eq.-9999.0d0) then
          phbetv=wbeta(4,is0)
        endif
      endif

      wlen=wtoe1/freqlow/1.0d9

      if (iundulator.ne.2) then
        sigrp=sqrt(wlen/(sourceeo(1,1,1)-sourceao(1,1,1)))
      else
        sigrp=sqrt(wlen/(dble(kampli)*phrperl))
      endif
      sigr=wlen/twopi1/sigrp
      !call util_break
      if (phwid.eq.-9999.0d0) then
        phwid=10.0d0*sqrt(sigr**2+
     &    ((phcenx-sourcen(1,1,1))*sigrp)**2)
      endif
      if (phhig.eq.-9999.0d0) then
        phhig=10.0d0*sqrt(sigr**2+
     &    ((phcenx-sourcen(1,1,1))*sigrp)**2)
      endif

      nphasez=(nphasez/2)*2+1
      nphasey=(nphasey/2)*2+1

      IF (NPHASEZ.GT.1) THEN
        DMASHZ=PHWID/(NPHASEZ-1)
      ELSE
        PHWID=0.D0
        DMASHZ=0.D0
      ENDIF
      IF (NPHASEY.GT.1) THEN
        DMASHY=PHHIG/(NPHASEY-1)
      ELSE
        PHHIG=0.D0
        DMASHY=0.D0
      ENDIF

      DA=PINW/max(1,(MOBSVZ-1))*PINH/max(1,(MOBSVY-1))

      if (phceny.eq.-9999.) phceny=ystart+vyin/vxin*(phcenx-xstart)
      if (phcenz.eq.-9999.) phcenz=zstart+vzin/vxin*(phcenx-xstart)

      if (iphfold.ne.0) then

        if (phbeth.le.0.0d0.or.phbetv.le.0.0d0) then
          write(lungfo,*)' '
          write(lungfo,*)
     &      '*** PHBETH or PHBETV lower or equal zero, WAVE aborted ***'
          write(lungfo,*)' '
          write(6,*)' '
          write(6,*)
     &      '*** PHBETH or PHBETV lower or equal zero, WAVE aborted ***'
          write(6,*)' '
        endif

        phgsigz=eps0h*phbeth
        phgsigy=eps0v*phbetv

        if (phgsigz.gt.0.0d0) then
          phgsigz=sqrt(phgsigz)
        else
          phgsigz=0.0d0
        endif

        if (phgsigy.gt.0.0d0) then
          phgsigy=sqrt(phgsigy)
        else
          phgsigy=0.0d0
        endif

      endif

      if (mphasez.eq.-9999) then
        if (dmashz.gt.0.0d0) then
          mphasez=nphasez+(dgsigz(1)*phgsigz/dmashz+1)*4
        else
          mphasez=0
        endif
      endif

      if (mphasey.eq.-9999) then
        if (dmashy.gt.0.0d0) then
          mphasey=nphasey+(dgsigy(1)*phgsigy/dmashy+1)*4
        else
          mphasey=0
        endif
      endif

      if (mphasez.lt.nphasez) mphasez=nphasez
      if (mphasey.lt.nphasey) mphasey=nphasey

      mphasez=(mphasez/2)*2+1
      mphasey=(mphasey/2)*2+1

      ALLOCATE(WSUM(NSOURCE*NFREQ))
      ALLOCATE(PHMEANZ(NSOURCE*NFREQ))
      ALLOCATE(PHMEANY(NSOURCE*NFREQ))
      ALLOCATE(PHSIGZ(NSOURCE*NFREQ))
      ALLOCATE(PHSIGY(NSOURCE*NFREQ))
      ALLOCATE(PHSHIFT(NOBSV))
      ALLOCATE(AMPLI(6,mphasez,mphasey,NFREQ))
      ALLOCATE(phspec3(mphasez,mphasey,nfreq))
      ALLOCATE(phspec3f(mphasez,mphasey,nfreq))
      ALLOCATE(phspec3fy(mphasez,mphasey,nfreq))
      ALLOCATE(zphw(mphasez))
      ALLOCATE(specwz(mphasez))
      ALLOCATE(specfwz(mphasez))
      ALLOCATE(yphw(mphasey))
      ALLOCATE(specwy(mphasey))
      ALLOCATE(specfwy(mphasey))
      allocate(phcoef(max(mphasez,mphasey)))
      allocate(phws1(max(mphasez,mphasey)))
      allocate(phws2(max(mphasez,mphasey)))
      allocate(phws3(max(mphasez,mphasey)))
      allocate(phws4(max(mphasez,mphasey)))
      ALLOCATE(EXPOM(NOBSV*NFREQ))
      ALLOCATE(DEXPOM(NOBSV))
      allocate(freq_omp(nfreq))

      ampli=(0.0d0,0.0d0)
      expom=(0.0d0,0.0d0)
      dexpom=(0.0d0,0.0d0)
      phshift=(0.0d0,0.0d0)
      specwz=0.0d0
      specwy=0.0d0
      phspec3=0.0d0
      phspec3f=0.0d0
      phspec3fy=0.0d0

      if (mhbookp.eq.0 .and. iroottrees.ge.0) then
        CALL hcdirm(OLDDIR,'R')
        CALL hropenm(LUNPH,'PHASE',FILEPH,'N',1024,ISTAT)
        CALL hcdirm(chphase(1:lenchphase),' ')
        IF (ISTAT.NE.0) THEN
          WRITE(6,*)'*** ERROR IN hropenm (PHASE_OMP) ***'
          WRITE(LUNGFO,*)'*** ERROR IN hropenm (PHASE_OMP) ***'
          STOP
        ENDIF
        CALL MHROUT(IDCODE,ICYCLE,' ')
      endif

      PHGEOSUM=0.0d0
      PHGEOSEL=0.0d0
      PHBEAM=0.0d0
      WSUM=0.0d0
      PHMEANZ=0.0d0
      PHMEANY=0.0d0
      PHSIGZ=0.0d0
      PHSIGY=0.0d0

      CALL hbookm(NIDPHASE,'PHASE',NTUP_P,chphase(1:lenchphase),
     &  mphasez*mphasey*nfreq,CHTAGS)

      if (nsource.gt.1) then
      CALL hbookm(NIDPHASE+1,'Field in pinhole (first source)',NTUP_P,chphase(1:lenchphase),
     &  nobsv*nfreq,CHTAGS)
      else
        CALL hbookm(NIDPHASE+1,'Field in pinhole',NTUP_P,chphase(1:lenchphase),
     &    nobsv*nfreq,CHTAGS)
      endif

      if (user(1).eq.1.0d0) then
        rea(1:2)=(0.0d0,0.0d0)
        rea(3)=dcmplx(reaima(3,1,icbrill),reaima(3,2,icbrill))
        reaima=0.0d0
        reaima(3,1,icbrill)=dreal(rea(3))
        reaima(3,2,icbrill)=dimag(rea(3))
      endif

      if (abs(phgshift).eq.9999.0d0) then
        do ifrq=1,nfreq
          iobfr=icbrill+nobsv*(ifrq-1)
          rea(1:2)=(0.0d0,0.0d0)
          rea(3)=dcmplx(reaima(3,1,iobfr),reaima(3,2,iobfr))
          expsh=rea(3)/abs(rea(3))
          if (phgshift.eq.-9999.0d0) expsh=expsh*cdexp(dcmplx(0.0d0,-pi1/2.0d0))
          DO iobs=1,nobsv
            iobfr=iobs+nobsv*(ifrq-1)
            rea=dcmplx(reaima(1:3,1,iobfr),reaima(1:3,2,iobfr))/expsh
            reaima(1:3,1,iobfr)=dreal(rea)
            reaima(1:3,2,iobfr)=dimag(rea)
          enddo
        enddo
      else if (phgshift.ne.0.0d0) then
        expsh=cdexp(dcmplx(0.0d0,phgshift))
        do ifrq=1,nfreq
          DO iobs=1,nobsv
            iobfr=iobs+nobsv*(ifrq-1)
            rea=dcmplx(reaima(1:3,1,iobfr),reaima(1:3,2,iobfr))/expsh
            reaima(1:3,1,iobfr)=dreal(rea)
            reaima(1:3,2,iobfr)=dimag(rea)
          enddo
        enddo
      endif

      isour=1
      smax=0.0d0
      do ifrq=1,nfreq
        DO iobs=1,nobsv
          iobfr=iobs+nobsv*(ifrq-1)
          if (spec(iobfr).gt.smax) then
            smax=spec(iobfr)
            reanor=
     &        reaima(1,1,iobfr)**2+reaima(1,2,iobfr)**2+
     &        reaima(2,1,iobfr)**2+reaima(2,2,iobfr)**2+
     &        reaima(3,1,iobfr)**2+reaima(3,2,iobfr)**2
          endif
        enddo
      enddo

      reanor=sqrt(smax/reanor)

      do ifrq=1,nfreq
        DO iobs=1,nobsv
          iobfr=iobs+nobsv*(ifrq-1)
          iphy=(iobs-1)/nobsvz+1
          iphz=iobs-(iphy-1)*nobsvz
c          if (idebug.ne.0.and.iphy.eq.nobsvy/2+1.and.ifrq.eq.2) then
c            write(99,*)iundulator,iphz,obsv(3,iobs),reaima(3,1,iobfr)*reanor
c          endif
          TUP(1:3)=obsv(1:3,iobs)
          if(abs(tup(2)).lt.1.0d-15) tup(2)=0.0d0
          if(abs(tup(3)).lt.1.0d-15) tup(3)=0.0d0
          TUP(4)=FREQ(ifrq)
          TUP(5)=ifrq
          TUP(6)=IPHY
          TUP(7)=IPHZ
          ef(1:3)=reanor*reaima(1:3,1,iobfr)
c          bf(1:3)=reanor*reaima(6:8,1,iobfr)
          efc(1:3)=dcmplx(reaima(1:3,1,iobfr),reaima(1:3,2,iobfr))
          bfc(1:3)=dcmplx(reaima(6:8,1,iobfr),reaima(6:8,2,iobfr))
          TUP(8)=ef(1)
          TUP(9)=reanor*reaima(1,2,iobfr)
          TUP(10)=ef(2)
          TUP(11)=reanor*reaima(2,2,iobfr)
          TUP(12)=ef(3)
          TUP(13)=reanor*reaima(3,2,iobfr)
          TUP(14)=spec(iobfr)
          if (ifold.ne.0) then
            TUP(15)=specf(iobfr)
          else
            tup(15)=0.0d0
          endif
          TUP(16)=reanor*reaima(6,1,iobfr)
          TUP(17)=reanor*reaima(6,2,iobfr)
          TUP(18)=reanor*reaima(7,1,iobfr)
          TUP(19)=reanor*reaima(7,2,iobfr)
          TUP(20)=reanor*reaima(8,1,iobfr)
          TUP(21)=reanor*reaima(8,2,iobfr)
c          rnx=ef(2)*bf(3)-ef(3)*bf(2)
c          rny=ef(3)*bf(1)-ef(1)*bf(3)
c          rnz=ef(1)*bf(2)-ef(2)*bf(1)
          rnx=real(efc(2)*conjg(bfc(3))-efc(3)*conjg(bfc(2)))
          rny=real(efc(3)*conjg(bfc(1))-efc(1)*conjg(bfc(3)))
          rnz=real(efc(1)*conjg(bfc(2))-efc(2)*conjg(bfc(1)))
          rn=sqrt(rnx**2+rny**2+rnz**2)
          if (rn.eq.0.0d0) rn=1.0d0
          tup(22)=rnx/rn
          tup(23)=rny/rn
          tup(24)=rnz/rn
          CALL hfm(NIDPHASE+1,TUP)
        ENDDO !iobs
      enddo !nfreq

      XPH=PHCENX
      XOBS=PINCEN(1)
      DX=XOBS-XPH
      DX2=DX*DX

      IF (DX.EQ.0.0D0) THEN
        WRITE(LUNGFO,*)'*** ERROR IN PHASE_OMP: PHCENX=PINCEN(1)  ***'
        WRITE(LUNGFO,*)'CHECK INPUT FILE'
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(6,*)'*** ERROR IN PHASE_OMP: PHCENX=PINCEN(1)  ***'
        WRITE(6,*)'CHECK INPUT FILE'
        WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
        STOP
      ENDIF

      PHLOWZ=PHCENZ-PHWID/2.D0
      PHLOWY=PHCENY-PHHIG/2.D0
      YPH=PHLOWY-DMASHY

      OMC=FREQ(1)/(HBAREV1*CLIGHT1)
      IF (ifreq2P.GT.2) THEN
        DOMC=(FREQ(2)-FREQ(1))/(HBAREV1*CLIGHT1)
      ELSE
        DOMC=OMC
      ENDIF !(ifrq2P.GT.2)

      mphasey_omp=mphasey
      mphasez_omp=mphasez
      nphasey_omp=nphasey
      nphasez_omp=nphasez
      iphfold_omp=iphfold
      nphelem_omp=nphelem
      ihsel_omp=ihsel
      nfreq_omp=nfreq
      iphase_omp=iphase
      iphfold_omp=iphfold

      phaperzm_omp=phaperzm
      phaperzp_omp=phaperzp
      phaperzpm_omp=phaperzpm
      phaperzpp_omp=phaperzpp
      phaperym_omp=phaperym
      phaperyp_omp=phaperyp
      phaperypm_omp=phaperypm
      phaperypp_omp=phaperypp
      phelem_omp=phelem
      freq_omp(1:nfreq)=freq(1:nfreq)

      iy1=(mphasey_omp-nphasey_omp)/2
      iz1=(mphasez_omp-nphasez_omp)/2

      dgsigy_omp=dgsigy(1)
      dgsigz_omp=dgsigz(1)

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMp& SHARED(mphasey_omp,mphasez_omp,nphasey_omp,nphasez_omp,iphfold_omp)
!$OMp& SHARED(nphelem_omp,ihsel_omp,nfreq_omp,freq_omp,wtoe1,iphase_omp)
!$OMP& SHARED(phaperzm_omp,phaperzp_omp,phaperzpm_omp,phaperzpp_omp)
!$OMP& SHARED(phaperym_omp,phaperyp_omp,phaperypm_omp,phaperypp_omp)
!$OMP& SHARED(phelem_omp,dmashy,dmashz,obsv,phlowy,phlowz,dx,dx2)
!$OMP& SHARED(ampli,reaima,domc,omc,da,obsvz,obsvy,yphw,zphw)

!$OMP& FIRSTPRIVATE(phshift,expom,dexpom,nobsvz,nobsvy,nobsv,iy1,iz1)
!$OMP& FIRSTPRIVATE(yai1y,yai2y,yai3y,yai4y,yai5y,yai6y)
!$OMP& FIRSTPRIVATE(yar1y,yar2y,yar3y,yar4y,yar5y,yar6y)
!$OMP& FIRSTPRIVATE(yai1,yai2,yai3,yai4,yai5,yai6)
!$OMP& FIRSTPRIVATE(yar1,yar2,yar3,yar4,yar5,yar6)
!$OMP& FIRSTPRIVATE(XA,xay,RESULT,RESULTY)

      allocate(
     &  x_th(max(nobsvy,nobsvz)),
     &  wobsv1_th(max(nobsvy,nobsvz)),wobsv2_th(max(nobsvy,nobsvz)),
     &  wobsv3_th(max(nobsvy,nobsvz)),wobsv4_th(max(nobsvy,nobsvz)),
     &  wobsv5_th(max(nobsvy,nobsvz)),wobsv6_th(max(nobsvy,nobsvz)),
     &  wobsv7_th(max(nobsvy,nobsvz)))

!$OMP DO

      DO iy=1,nphasey_omp

        iphy=iy+iy1

        ith=OMP_GET_THREAD_NUM()+1
        yph=phlowy+dble(iy-1)*dmashy

        DO iz=1,nphasez_omp

          iphz=iz+iz1

          zph=phlowz+dble(iz-1)*dmashz

          DO IOBS=1,NOBSV

            XOBS=OBSV(1,IOBS)
            YOBS=OBSV(2,IOBS)
            ZOBS=OBSV(3,IOBS)

            DY=YOBS-YPH
            DZ=ZOBS-ZPH
            DZY2=DZ*DZ+DY*DY

C     TO MAKE SURE THAT TAYLOR-EXPANSION IS VALID

            IF (DZY2.le.0.01D0*DX2) THEN

              EPS(1)=DZY2/DX2
              DO IEPS=2,6
                EPS(IEPS)=EPS(IEPS-1)*EPS(1)
              ENDDO !IEPS

c      TAYLOR-EXPANSION DONE WITH REDUCE
c     IN "WTAY1.RED";
c     on rounded;
c     on numval;
c     precision 13;
c     F:=SQRT(1+EPS);
c     DR:=TAY1(F,EPS,6);
c     ON FORT;
c     OUT "RED.FOR";
c     DR;
c     SHUT "RED.FOR";
C ans is actually reduce by 1.0 to avoid large overall phase

              ans=-0.0205078125D0*eps(6)+0.02734375D0*eps(5)
     &          -0.0390625D0*eps(4)+
     &          0.0625D0*eps(3)-0.125D0*eps(2)+0.5D0*eps(1)

              DR=DABS(DX*(ANS+1.D0))
              DRRED=-DABS(DX*ANS)

            else
c              WRITE(LUNGFO,*)'*** ERROR IN PHASE_OMP: DZY2.GT.0.01D0*DX2  ***'
c              WRITE(LUNGFO,*)'CHECK INPUT FILE AND INCREASE PINCEN(1)'
c              WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
c              WRITE(6,*)'*** ERROR IN PHASE_OMP: PHCENX=PINCEN(1)  ***'
c              WRITE(6,*)'CHECK INPUT FILE AND INCREASE PINCEN(1)'
c              WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
c              STOP
              dr=sqrt(1.0d0+dzy2/dx2)*abs(dx)
              drred=-dabs(dr-abs(dx))
            ENDIF

            IF (DR.NE.0.0d0) THEN
              EXPOM(IOBS)=CDEXP(DCMPLX(0.0d0,DRRED*OMC))/DR
            ELSE
              EXPOM(IOBS)=1.0D0
            ENDIF

            DEXPOM(IOBS)=CDEXP(DCMPLX(0.0d0,DRRED*DOMC))
c            print*,ith,iobs,expom(iobs)
c+seq,dum2.
          ENDDO   !NOBS

          DO ifrq=1,NFREQ_omp

            RLAMBDA1=FREQ_omp(ifrq)/WTOE1*1.0D9   !1/lambda[m]=1/(wtoe1/freq*1.e-9)

            IF (IPHASE_omp.GT.0) THEN

              DO IOBS=1,NOBSV

                IOBFR=IOBS+NOBSV*(ifrq-1)

                IF (ifrq.EQ.1) THEN
                  PHSHIFT(IOBS)=EXPOM(IOBFR)
                ELSE
                  PHSHIFT(IOBS)=PHSHIFT(IOBS)*DEXPOM(IOBS)
                ENDIF   !(ifrq.EQ.1)
c                print*,iobfr,iobs,expom(iobfr),phshift(iobs)
c                stop
                IF (DX.GE.0) THEN

                  ampli(1,iphz,iphy,ifrq)=ampli(1,iphz,iphy,ifrq)+
     &              DCMPLX(REAIMA(1,1,IOBFR),REAIMA(1,2,IOBFR))
     &              *PHSHIFT(IOBS)
                  ampli(2,iphz,iphy,ifrq)=ampli(2,iphz,iphy,ifrq)+
     &              DCMPLX(REAIMA(2,1,IOBFR),REAIMA(2,2,IOBFR))
     &              *PHSHIFT(IOBS)
                  ampli(3,iphz,iphy,ifrq)=ampli(3,iphz,iphy,ifrq)+
     &              DCMPLX(REAIMA(3,1,IOBFR),REAIMA(3,2,IOBFR))
     &              *PHSHIFT(IOBS)

                  ampli(4,iphz,iphy,ifrq)=ampli(4,iphz,iphy,ifrq)+
     &              DCMPLX(reaima(6,1,IOBFR),reaima(6,2,IOBFR))
     &              *PHSHIFT(IOBS)
                  ampli(5,iphz,iphy,ifrq)=ampli(5,iphz,iphy,ifrq)+
     &              DCMPLX(reaima(7,1,IOBFR),reaima(7,2,IOBFR))
     &              *PHSHIFT(IOBS)
                  ampli(6,iphz,iphy,ifrq)=ampli(6,iphz,iphy,ifrq)+
     &              DCMPLX(reaima(8,1,IOBFR),reaima(8,2,IOBFR))
     &              *PHSHIFT(IOBS)

                ELSE

                  ampli(1,iphz,iphy,ifrq)=ampli(1,iphz,iphy,ifrq)+
     &              DCMPLX(REAIMA(1,1,IOBFR),-REAIMA(1,2,IOBFR))
     &              *PHSHIFT(IOBS)
                  ampli(2,iphz,iphy,ifrq)=ampli(2,iphz,iphy,ifrq)+
     &              DCMPLX(REAIMA(2,1,IOBFR),-REAIMA(2,2,IOBFR))
     &              *PHSHIFT(IOBS)
                  ampli(3,iphz,iphy,ifrq)=ampli(3,iphz,iphy,ifrq)+
     &              DCMPLX(REAIMA(3,1,IOBFR),-REAIMA(3,2,IOBFR))
     &              *PHSHIFT(IOBS)

                  ampli(4,iphz,iphy,ifrq)=ampli(4,iphz,iphy,ifrq)+
     &              DCMPLX(reaima(6,1,IOBFR),-REAIMA(6,2,IOBFR))
     &              *PHSHIFT(IOBS)
                  ampli(5,iphz,iphy,ifrq)=ampli(5,iphz,iphy,ifrq)+
     &              DCMPLX(reaima(7,1,IOBFR),-reaima(7,2,IOBFR))
     &              *PHSHIFT(IOBS)
                  ampli(6,iphz,iphy,ifrq)=ampli(6,iphz,iphy,ifrq)+
     &              DCMPLX(reaima(8,1,IOBFR),-reaima(8,2,IOBFR))
     &              *PHSHIFT(IOBS)

                ENDIF !(DX.GE.0)

              ENDDO  !NOBSV

              ampli(1:6,iphz,iphy,ifrq)=ampli(1:6,iphz,iphy,ifrq)*DA*RLAMBDA1

            ELSE  !IPHASE_omp.GT.0

              DO IOBSY=1,NOBSVY
                DO IOBSZ=1,NOBSVZ

                  iobs=(iobsy-1)*nobsvz+iobsz

                  IF (ifrq.EQ.1) THEN
                    PHSHIFT(IOBS)=EXPOM(IOBS+NOBSV*(ifrq-1))
                  ELSE
                    PHSHIFT(IOBS)=PHSHIFT(IOBS)*DEXPOM(IOBS)
                  ENDIF !(ifrq.EQ.1)

                  IOBFR=IOBS+NOBSV*(ifrq-1)
c+seq,dummy.
                  IF (DX.GE.0) THEN

                    ampli(1,iphz,iphy,ifrq)=
     &                DCMPLX(REAIMA(1,1,IOBFR),REAIMA(1,2,IOBFR))
     &                *PHSHIFT(IOBS)
                    ampli(2,iphz,iphy,ifrq)=
     &                DCMPLX(REAIMA(2,1,IOBFR),REAIMA(2,2,IOBFR))
     &                *PHSHIFT(IOBS)
                    ampli(3,iphz,iphy,ifrq)=
     &                DCMPLX(REAIMA(3,1,IOBFR),REAIMA(3,2,IOBFR))
     &                *PHSHIFT(IOBS)

                    ampli(4,iphz,iphy,ifrq)=
     &                DCMPLX(reaima(6,1,IOBFR),-reaima(6,2,IOBFR))
     &                *PHSHIFT(IOBS)
                    ampli(5,iphz,iphy,ifrq)=
     &                DCMPLX(reaima(7,1,IOBFR),-reaima(7,2,IOBFR))
     &                *PHSHIFT(IOBS)
                    ampli(6,iphz,iphy,ifrq)=
     &                DCMPLX(reaima(8,1,IOBFR),-reaima(8,2,IOBFR))
     &                *PHSHIFT(IOBS)

                  ELSE

                    ampli(1,iphz,iphy,ifrq)=
     &                DCMPLX(REAIMA(1,1,IOBFR),-REAIMA(1,2,IOBFR))
     &                *PHSHIFT(IOBS)
                    ampli(2,iphz,iphy,ifrq)=
     &                DCMPLX(REAIMA(2,1,IOBFR),-REAIMA(2,2,IOBFR))
     &                *PHSHIFT(IOBS)
                    ampli(3,iphz,iphy,ifrq)=
     &                DCMPLX(REAIMA(3,1,IOBFR),-REAIMA(3,2,IOBFR))
     &                *PHSHIFT(IOBS)

                    ampli(4,iphz,iphy,ifrq)=
     &                DCMPLX(reaima(6,1,IOBFR),+reaima(6,2,IOBFR))
     &                *PHSHIFT(IOBS)
                    ampli(5,iphz,iphy,ifrq)=
     &                DCMPLX(reaima(7,1,IOBFR),+reaima(7,2,IOBFR))
     &                *PHSHIFT(IOBS)
                    ampli(6,iphz,iphy,ifrq)=
     &                DCMPLX(reaima(8,1,IOBFR),+reaima(8,2,IOBFR))
     &                *PHSHIFT(IOBS)

                  ENDIF !(DX.GE.0)

                  X_th(IOBSZ)=OBSVZ(IOBSZ)

                  YAR1(IOBSZ)=DREAL(ampli(1,iphz,iphy,ifrq))
                  YAR2(IOBSZ)=DREAL(ampli(2,iphz,iphy,ifrq))
                  YAR3(IOBSZ)=DREAL(ampli(3,iphz,iphy,ifrq))

                  YAI1(IOBSZ)=DIMAG(ampli(1,iphz,iphy,ifrq))
                  YAI2(IOBSZ)=DIMAG(ampli(2,iphz,iphy,ifrq))
                  YAI3(IOBSZ)=DIMAG(ampli(3,iphz,iphy,ifrq))

                  YAR4(IOBSZ)=DREAL(ampli(4,iphz,iphy,ifrq))
                  YAR5(IOBSZ)=DREAL(ampli(5,iphz,iphy,ifrq))
                  YAR6(IOBSZ)=DREAL(ampli(6,iphz,iphy,ifrq))

                  YAI4(IOBSZ)=DIMAG(ampli(4,iphz,iphy,ifrq))
                  YAI5(IOBSZ)=DIMAG(ampli(5,iphz,iphy,ifrq))
                  YAI6(IOBSZ)=DIMAG(ampli(6,iphz,iphy,ifrq))

                ENDDO   !NOBSVZ

                wobsv1_th(1:nobsvz)=yar1(1:nobsvz)
                CALL UTIL_SPLINE_integral_omp(nobsvz,RESULT(1,1))
                wobsv1_th(1:nobsvz)=yai1(1:nobsvz)
                CALL UTIL_SPLINE_integral_omp(nobsvz,RESULT(1,2))

                wobsv1_th(1:nobsvz)=yar2(1:nobsvz)
                CALL UTIL_SPLINE_integral_omp(nobsvz,RESULT(2,1))
                wobsv1_th(1:nobsvz)=yai2(1:nobsvz)
                CALL UTIL_SPLINE_integral_omp(nobsvz,RESULT(2,2))

                wobsv1_th(1:nobsvz)=yar3(1:nobsvz)
                CALL UTIL_SPLINE_integral_omp(nobsvz,RESULT(3,1))
                wobsv1_th(1:nobsvz)=yai3(1:nobsvz)
                CALL UTIL_SPLINE_integral_omp(nobsvz,RESULT(3,2))

                wobsv1_th(1:nobsvz)=yar4(1:nobsvz)
                CALL UTIL_SPLINE_integral_omp(nobsvz,RESULT(4,1))
                wobsv1_th(1:nobsvz)=yai4(1:nobsvz)
                CALL UTIL_SPLINE_integral_omp(nobsvz,RESULT(4,2))

                wobsv1_th(1:nobsvz)=yar5(1:nobsvz)
                CALL UTIL_SPLINE_integral_omp(nobsvz,RESULT(5,1))
                wobsv1_th(1:nobsvz)=yai5(1:nobsvz)
                CALL UTIL_SPLINE_integral_omp(nobsvz,RESULT(5,2))

                wobsv1_th(1:nobsvz)=yar6(1:nobsvz)
                CALL UTIL_SPLINE_integral_omp(nobsvz,RESULT(6,1))
                wobsv1_th(1:nobsvz)=yai6(1:nobsvz)
                CALL UTIL_SPLINE_integral_omp(nobsvz,RESULT(6,2))

                XAY(IOBSY)=OBSVY(IOBSY)

                YAR1Y(IOBSY)=RESULT(1,1)
                YAI1Y(IOBSY)=RESULT(1,2)
                YAR2Y(IOBSY)=RESULT(2,1)
                YAI2Y(IOBSY)=RESULT(2,2)
                YAR3Y(IOBSY)=RESULT(3,1)
                YAI3Y(IOBSY)=RESULT(3,2)
                YAR4Y(IOBSY)=RESULT(4,1)
                YAI4Y(IOBSY)=RESULT(4,2)
                YAR5Y(IOBSY)=RESULT(5,1)
                YAI5Y(IOBSY)=RESULT(5,2)
                YAR6Y(IOBSY)=RESULT(6,1)
                YAI6Y(IOBSY)=RESULT(6,2)

              ENDDO !NOBSVY

              x_th(1:nobsvy)=xay(1:nobsvy)

              wobsv1_th(1:nobsvy)=yar1y(1:nobsvy)
              CALL UTIL_SPLINE_integral_omp(NOBSVY,RESULTY(1,1))
              wobsv1_th(1:nobsvy)=yai1y(1:nobsvy)
              CALL UTIL_SPLINE_integral_omp(NOBSVY,RESULTY(1,2))

              wobsv1_th(1:nobsvy)=yar2y(1:nobsvy)
              CALL UTIL_SPLINE_integral_omp(NOBSVY,RESULTY(2,1))
              wobsv1_th(1:nobsvy)=yai2y(1:nobsvy)
              CALL UTIL_SPLINE_integral_omp(NOBSVY,RESULTY(2,2))

              wobsv1_th(1:nobsvy)=yar3y(1:nobsvy)
              CALL UTIL_SPLINE_integral_omp(NOBSVY,RESULTY(3,1))
              wobsv1_th(1:nobsvy)=yai3y(1:nobsvy)
              CALL UTIL_SPLINE_integral_omp(NOBSVY,RESULTY(3,2))

              wobsv1_th(1:nobsvy)=yar4y(1:nobsvy)
              CALL UTIL_SPLINE_integral_omp(NOBSVY,RESULTY(4,1))
              wobsv1_th(1:nobsvy)=yai4y(1:nobsvy)
              CALL UTIL_SPLINE_integral_omp(NOBSVY,RESULTY(4,2))

              wobsv1_th(1:nobsvy)=yar5y(1:nobsvy)
              CALL UTIL_SPLINE_integral_omp(NOBSVY,RESULTY(5,1))
              wobsv1_th(1:nobsvy)=yai5y(1:nobsvy)
              CALL UTIL_SPLINE_integral_omp(NOBSVY,RESULTY(5,2))

              wobsv1_th(1:nobsvy)=yar6y(1:nobsvy)
              CALL UTIL_SPLINE_integral_omp(NOBSVY,RESULTY(6,1))
              wobsv1_th(1:nobsvy)=yai6y(1:nobsvy)
              CALL UTIL_SPLINE_integral_omp(NOBSVY,RESULTY(6,2))

              ampli(1,iphz,iphy,ifrq)=DCMPLX(RESULTY(1,1),RESULTY(1,2))*RLAMBDA1
              ampli(2,iphz,iphy,ifrq)=DCMPLX(RESULTY(2,1),RESULTY(2,2))*RLAMBDA1
              ampli(3,iphz,iphy,ifrq)=DCMPLX(RESULTY(3,1),RESULTY(3,2))*RLAMBDA1

              ampli(4,iphz,iphy,ifrq)=DCMPLX(RESULTY(4,1),RESULTY(4,2))*RLAMBDA1
              ampli(5,iphz,iphy,ifrq)=DCMPLX(RESULTY(5,1),RESULTY(5,2))*RLAMBDA1
              ampli(6,iphz,iphy,ifrq)=DCMPLX(RESULTY(6,1),RESULTY(6,2))*RLAMBDA1

            ENDIF !IPHASE.GT.0

          ENDDO   !NFREQ

        ENDDO  !NPHASEZ

      ENDDO !NPHAZEY

!$OMP END DO

      deallocate(x_th,
     &  wobsv1_th,wobsv2_th,wobsv3_th,wobsv4_th,wobsv5_th,
     &  wobsv6_th,wobsv7_th)

!$OMP END PARALLEL

      do iphz=1,mphasez
        zphw(iphz)=-(mphasez-1)*dmashz/2.0+(iphz-1)*dmashz
      enddo

      do iphy=1,mphasey
        yphw(iphy)=-(mphasey-1)*dmashy/2.0+(iphy-1)*dmashy
      enddo

      !print*,iz1,iy1,dgsigz_omp,dgsigy_omp

      smax=-1.0d30

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMp& SHARED(mphasey_omp,mphasez_omp,nphasey_omp,nphasez_omp,iphfold_omp)
!$OMp& SHARED(nphelem_omp,ihsel_omp,nfreq_omp,freq_omp,wtoe1,iphase_omp)
!$OMP& SHARED(phaperzm_omp,phaperzp_omp,phaperzpm_omp,phaperzpp_omp)
!$OMP& SHARED(phaperym_omp,phaperyp_omp,phaperypm_omp,phaperypp_omp)
!$OMP& SHARED(phelem_omp,dmashy,dmashz,obsv,phlowy,phlowz,dx,dx2)
!$OMP& SHARED(ampli,reaima,phshift,domc,omc,da,obsvz,obsvy,yphw,zphw)
!$OMP& SHARED(smax,phspec3,sfmax,phspec3fy,phspec3f,dgsigz_omp,dgsigy_omp,phgsigz,phgsigy)

!$OMP& FIRSTPRIVATE(nobsvz,nobsvy,nobsv,nidphase,xph,iy1,iz1)
!$OMP& FIRSTPRIVATE(specwy,specfwy,specwz,specfwz)

      allocate(
     &  x_th(max(mphasez_omp,mphasey_omp)),
     &  wobsv1_th(max(mphasez_omp,mphasey_omp)),wobsv2_th(max(mphasez_omp,mphasey_omp)),
     &  wobsv3_th(max(mphasez_omp,mphasey_omp)),wobsv4_th(max(mphasez_omp,mphasey_omp)),
     &  wobsv5_th(max(mphasez_omp,mphasey_omp)),wobsv6_th(max(mphasez_omp,mphasey_omp)),
     &  wobsv7_th(max(mphasez_omp,mphasey_omp)))

!$OMP DO

      do ifrq=1,nfreq_omp

        DO iy=1,nphasey_omp
          iphy=iy+iy1
          DO iz=1,nphasez_omp
            iphz=iz+iz1
            phspec3(iphz,iphy,ifrq)=
     &        DREAL(ampli(1,iphz,iphy,ifrq))*DREAL(ampli(1,iphz,iphy,ifrq))+
     &        DIMAG(ampli(1,iphz,iphy,ifrq))*DIMAG(ampli(1,iphz,iphy,ifrq))+
     &        DREAL(ampli(2,iphz,iphy,ifrq))*DREAL(ampli(2,iphz,iphy,ifrq))+
     &        DIMAG(ampli(2,iphz,iphy,ifrq))*DIMAG(ampli(2,iphz,iphy,ifrq))+
     &        DREAL(ampli(3,iphz,iphy,ifrq))*DREAL(ampli(3,iphz,iphy,ifrq))+
     &        DIMAG(ampli(3,iphz,iphy,ifrq))*DIMAG(ampli(3,iphz,iphy,ifrq))
            if (phspec3(iphz,iphy,ifrq).gt.smax) smax=phspec3(iphz,iphy,ifrq)
          ENDDO   !NPHASez_omp
        ENDDO  !NPHAZey_omp

        if (iphfold_omp.ne.0) then

          do iphy=(mphasey_omp-nphasey_omp)/2+1,(mphasey_omp-nphasey_omp)/2+nphasey_omp

            specwz=phspec3(1:mphasez_omp,iphy,ifrq)

            if (dgsigz_omp.gt.0.0d0.and.phgsigz.gt.0.0d0) then
              if (iphfold_omp.lt.0) then
                x_th(1:mphasez_omp)=zphw(1:mphasez_omp)
                wobsv1_th(1:mphasez_omp)=specwz(1:mphasez_omp)
                call util_fold_function_gauss_lin_omp(
     &            mphasez_omp,phgsigz,dgsigz_omp)
                specfwz(1:mphasez_omp)=wobsv2_th(1:mphasez_omp)
              else
                x_th(1:mphasez_omp)=zphw(1:mphasez_omp)
                wobsv1_th(1:mphasez_omp)=specwz(1:mphasez_omp)
                call util_fold_function_gauss_omp(
     &            mphasez_omp,phgsigz,dgsigz_omp)
                specfwz(1:mphasez_omp)=wobsv2_th(1:mphasez_omp)
              endif
            else
              specfwz=specwz
            endif

            phspec3fy(1:mphasez_omp,iphy,ifrq)=specfwz(1:mphasez_omp)

          enddo !iphy

          do iphz=1,mphasez_omp
            specwy=phspec3fy(iphz,1:mphasey_omp,ifrq)
            if (dgsigy_omp.gt.0.0d0.and.phgsigy.gt.0.0d0) then
              if (iphfold_omp.lt.0) then
c               call util_fold_function_gauss_lin(
c     &            mphasey_omp,yphw,specwy,phgsigy,dgsigy_omp,specfwy,phws1,phws2)
                x_th(1:mphasey_omp)=yphw(1:mphasey_omp)
                wobsv1_th(1:mphasey_omp)=specwy(1:mphasey_omp)
                call util_fold_function_gauss_lin_omp(
     &            mphasey_omp,phgsigy,dgsigy_omp)
                specfwy(1:mphasey_omp)=wobsv2_th(1:mphasey_omp)
              else
c                call util_fold_function_gauss(
c     &            mphasey_omp,yphw,specwy,phgsigy,dgsigy_omp,specfwy,phcoef,
c     &            phws1,phws2,phws3,phws4)
                x_th(1:mphasey_omp)=yphw(1:mphasey_omp)
                wobsv1_th(1:mphasey_omp)=specwy(1:mphasey_omp)
                call util_fold_function_gauss_omp(
     &            mphasey_omp,phgsigy,dgsigy_omp)
                specfwy(1:mphasey_omp)=wobsv2_th(1:mphasey_omp)
              endif
            else
              specfwy=specwy
            endif
            phspec3f(iphz,1:mphasey_omp,ifrq)=specfwy
          enddo !iphz

        endif !iphfold_omp

      enddo !nfreq

!$OMP END DO

      deallocate(x_th,
     &  wobsv1_th,wobsv2_th,wobsv3_th,wobsv4_th,wobsv5_th,
     &  wobsv6_th,wobsv7_th)

!$OMP END PARALLEL

      sfmax=-1.0d30

      do ifrq=1,nfreq
        DO iphy=1,mphasey_omp
          DO IPHZ=1,mphasez_omp
            TUP(1)=XPH
            TUP(2)=YPHw(iphy)
            TUP(3)=ZPHw(iphz)
            if(abs(tup(2)).lt.1.0d-15) tup(2)=0.0d0
            if(abs(tup(3)).lt.1.0d-15) tup(3)=0.0d0
            TUP(4)=FREQ_omp(ifrq)
            TUP(5)=ifrq
            TUP(6)=IPHY
            TUP(7)=IPHZ
            ef(1:3)=DREAL(reanor*ampli(1:3,iphz,iphy,ifrq))
c            bf(1:3)=DREAL(reanor*ampli(4:6,iphz,iphy,ifrq))
            efc(1:3)=ampli(1:3,iphz,iphy,ifrq)
            bfc(1:3)=ampli(4:6,iphz,iphy,ifrq)
            TUP(8)=ef(1)
            TUP(9)=DIMAG(reanor*ampli(1,iphz,iphy,ifrq))
            TUP(10)=ef(2)
            TUP(11)=DIMAG(reanor*ampli(2,iphz,iphy,ifrq))
            TUP(12)=ef(3)
            TUP(13)=DIMAG(reanor*ampli(3,iphz,iphy,ifrq))
            TUP(14)=phspec3(iphz,iphy,ifrq)*reanor**2
            TUP(15)=phspec3f(iphz,iphy,ifrq)*reanor**2
            TUP(16)=DREAL(reanor*ampli(4,iphz,iphy,ifrq))
            TUP(17)=DIMAG(reanor*ampli(4,iphz,iphy,ifrq))
            TUP(18)=DREAL(reanor*ampli(5,iphz,iphy,ifrq))
            TUP(19)=DIMAG(reanor*ampli(5,iphz,iphy,ifrq))
            TUP(20)=DREAL(reanor*ampli(6,iphz,iphy,ifrq))
            TUP(21)=DIMAG(reanor*ampli(6,iphz,iphy,ifrq))
c            rnx=ef(2)*bf(3)-ef(3)*bf(2)
c            rny=ef(3)*bf(1)-ef(1)*bf(3)
c            rnz=ef(1)*bf(2)-ef(2)*bf(1)
            rnx=real(efc(2)*conjg(bfc(3))-efc(3)*conjg(bfc(2)))
            rny=real(efc(3)*conjg(bfc(1))-efc(1)*conjg(bfc(3)))
            rnz=real(efc(1)*conjg(bfc(2))-efc(2)*conjg(bfc(1)))
            rn=sqrt(rnx**2+rny**2+rnz**2)
            if (rn.eq.0.0d0) rn=1.0d0
            tup(22)=rnx/rn
            tup(23)=rny/rn
            tup(24)=rnz/rn
            if (phspec3f(iphz,iphy,ifrq).gt.sfmax) sfmax=phspec3f(iphz,iphy,ifrq)
            CALL hfm(NIDPHASE,TUP)
          ENDDO   !mPHASez_omp
        ENDDO  !MPHASey_omp
      enddo !nfreq

      if (iphfold.ne.0.and.sfmax.ge.smax) then
        write(lungfo,*)' '
        write(lungfo,*)
     &    '*** Warning in PHASE: Max. of folded intensity higher than unfolded one!'
        write(lungfo,*)'   Be careful with results, they are probably wrong'
        write(lungfo,*)'   Try different parameters of tiny beam current to investigate the problem.'
        write(lungfo,*)'   Max. of raw and folded intensities:'
        write(lungfo,*)'   ',smax,sfmax
        write(lungfo,*)' '
        write(6,*)' '
        write(6,*)
     &    '*** Warning in PHASE: Max. of folded intensity higher than unfolded one!'
        write(6,*)'   Be careful with results, they are probably wrong'
        write(6,*)'   Try different parameters of tiny beam current to investigate the problem.'
        write(6,*)'   Max. of raw and folded intensities:'
        write(6,*)'   ',smax,sfmax
        write(6,*)' '
      endif

      CALL MHROUT(NIDPHASE,ICYCLE,' ')
      CALL hdeletm(NIDPHASE)

      IF (mPHASEZ.GT.1.AND.mPHASEY.GT.1) THEN

        call hbook2m(NIDPHASE-1,'PHASE',
     &    mPHASEZ,
     &    SNGL(PHCENZ-(mphasez-1)*dmashz/2.-PHWID/(NPHASEZ-1)/2.),
     &    SNGL(PHCENZ+(mphasez-1)*dmashz/2.+PHWID/(NPHASEZ-1)/2.),
     &    mPHASEY,
     &    SNGL(PHCENY-(mphasey-1)*dmashy/2.-PHHIG/(NPHASEY-1)/2.),
     &    SNGL(PHCENY+(mphasey-1)*dmashy/2.+PHHIG/(NPHASEY-1)/2.),
     &    0.0)
        CALL MHROUT(NIDPHASE-1,ICYCLE,' ')
        CALL hdeletm(NIDPHASE-1)

        call hbook1m(NIDPHASE-2,'PHASE (HORIZONTAL CUT)',
     &    mPHASEZ,
     &    SNGL(PHCENZ-(mphasez-1)*dmashz/2.-PHWID/(NPHASEZ-1)/2.),
     &    SNGL(PHCENZ+(mphasez-1)*dmashz/2.+PHWID/(NPHASEZ-1)/2.),
     &    0.0)
        CALL MHROUT(NIDPHASE-2,ICYCLE,' ')
        CALL hdeletm(NIDPHASE-2)

        call hbook1m(NIDPHASE-3,'PHASE (VERTICAL CUT)',
     &    mPHASEY,
     &    SNGL(PHCENY-(mphasey-1)*dmashy/2.-PHHIG/(NPHASEY-1)/2.),
     &    SNGL(PHCENY+(mphasey-1)*dmashy/2.+PHHIG/(NPHASEY-1)/2.),
     &    0.0)
        CALL MHROUT(NIDPHASE-3,ICYCLE,' ')
        CALL hdeletm(NIDPHASE-3)

      ELSE IF (NPHASEZ.GT.1) THEN

        call hbook1m(NIDPHASE-2,'PHASE (HORIZONTAL CUT)',
     &    mPHASEZ,
     &    SNGL(PHCENZ-(mphasez-1)*dmashz/2.-PHWID/(NPHASEZ-1)/2.),
     &    SNGL(PHCENZ+(mphasez-1)*dmashz/2.+PHWID/(NPHASEZ-1)/2.),
     &    0.0)
        CALL MHROUT(NIDPHASE-2,ICYCLE,' ')
        CALL hdeletm(NIDPHASE-2)

      ELSE IF (NPHASEY.GT.1) THEN

        call hbook1m(NIDPHASE-3,'PHASE (VERTICAL CUT)',
     &    mPHASEY,
     &    SNGL(PHCENY-(mphasey-1)*dmashy/2.-PHHIG/(NPHASEY-1)/2.),
     &    SNGL(PHCENY+(mphasey-1)*dmashy/2.+PHHIG/(NPHASEY-1)/2.),
     &    0.0)
        CALL MHROUT(NIDPHASE-3,ICYCLE,' ')
        CALL hdeletm(NIDPHASE-3)

      ENDIF !(NPHASEZ.GT.1.AND.NPHASEY.GT.1) THEN

      if (mhbookp.eq.0.and.iroottrees.ge.0) then
        CALL hrendm('PHASE')
        CLOSE(LUNPH)
      endif

      if (mhbookp.eq.0.and.iroottrees.ge.0) then
        CALL hropenm(LUNPH,'PHASE','GEO_'//FILEPH,'N',1024,ISTAT)
        CALL hcdirm(chphase(1:lenchphase),' ')
        IF (ISTAT.NE.0) THEN
          WRITE(6,*)'*** ERROR IN hropenm (PHASE_OMP) ***'
          WRITE(LUNGFO,*)'*** ERROR IN hropenm (PHASE_OMP) ***'
          STOP
        ENDIF
        CALL MHROUT(IDCODE,ICYCLE,' ')
      endif

      IF (ABS(IPHASE).GT.1) THEN

C PHASE SPACE FORM GEOMETRICAL OPTIC. ONLY CORRECT IF DIFFRACTION IS
C NEGLIGIBLE (TO BE CHECKED BY TRANSFORMED PHASE)
        CALL hbookm(NIDGEO,'PHASE SPACE DIST. (GEO. OPTIC)'
     &    ,NGEO_P,chphase(1:lenchphase),nphelem*nsource*nobsv,CHGEO)

        NIDGEO1=NIDGEO+1
        CALL hbookm(NIDGEO1,'SEL. PHASE SPACE DIST. (GEO. OPTIC)'
     &    ,NGEO_P,chphase(1:lenchphase),nphelem*nsource*nobsv*nfreq,CHGEO)

        NIDGEO2=NIDGEO+2
        CALL hbookm(NIDGEO2,
     &    'SELECTED PHASE SPACE DIST. AT END OF BEAMLINE'
     &    ,NBEAM_P,chphase(1:lenchphase),nphelem*nsource*nobsv*nfreq,CHBEAM)

        CALL hbookm(NIDGEO+3,'SOURCE SIZE (GEO. OPTIC THROUGH APERTURE)'
     &    ,NSIZE_P,chphase(1:lenchphase),nphelem*nsource*nobsv*nfreq,CHSIZE)
C--- GET FOCUSSING

        DO J=1,4
          DO I=1,4
            TOTMAT(I,J)=0.0D0
            DUMMAT(I,J)=0.0D0
          ENDDO
        ENDDO

        TOTMAT(1,1)=1.0D0
        TOTMAT(2,2)=1.0D0
        TOTMAT(3,3)=1.0D0
        TOTMAT(4,4)=1.0D0

        DO IELEM=1,NPHELEM
          DO J=1,4
            DO I=1,4
              DO K=1,4
                DUMMAT(I,J)=DUMMAT(I,J)+TOTMAT(I,K)*PHELEM(K,J,IELEM)
              ENDDO
            ENDDO
          ENDDO
          DO J=1,4
            DO I=1,4
              TOTMAT(I,J)=DUMMAT(I,J)
              DUMMAT(I,J)=0.0D0
            ENDDO
          ENDDO
        ENDDO   !IELEM

        ZBEAM=EPSBEAM
        TANPHIB=0.0
        YBEAM=EPSBEAM
        TANTHEB=0.0

        DO IELEM=1,NPHELEM

          BEAM(1)=ZBEAM
          BEAM(2)=TANPHIB
          BEAM(3)=YBEAM
          BEAM(4)=TANTHEB

          DO J=1,4
            DO I=1,4
              OPTMAT(I,J)=PHELEM(I,J,IELEM)
            ENDDO
          ENDDO

          IF (OPTMAT(2,1).GT.0.0) THEN
            WRITE(6,*)'*** WARNING IN PHASE: PHELEM(2,1,N).GT.0'
            WRITE(6,*)'PHELEM(2,1,N) IS -1/fx !'
            WRITE(6,*)'CHECK ELEMENT ',IELEM
            WRITE(LUNGFO,*)'*** WARNING IN PHASE: PHELEM(2,1,N).GT.0'
            WRITE(LUNGFO,*)'PHELEM(2,1,N) IS -1/fx !'
            WRITE(LUNGFO,*)'CHECK ELEMENT ',IELEM
          ENDIF

          IF (OPTMAT(4,3).GT.0.0) THEN
            WRITE(6,*)'*** WARNING IN PHASE: PHELEM(4,3,N).GT.0'
            WRITE(6,*)'PHELEM(4,3,N) IS -1/fx !'
            WRITE(6,*)'CHECK ELEMENT ',IELEM
            WRITE(LUNGFO,*)'*** WARNING IN PHASE: PHELEM(2,1,N).GT.0'
            WRITE(LUNGFO,*)'PHELEM(4,3,N) IS -1/fy !'
            WRITE(LUNGFO,*)'CHECK ELEMENT ',IELEM
          ENDIF

          ZBEAM=
     &      OPTMAT(1,1)*BEAM(1)
     &      +OPTMAT(1,2)*BEAM(2)
     &      +OPTMAT(1,3)*BEAM(3)
     &      +OPTMAT(1,4)*BEAM(4)

          TANPHIB=
     &      OPTMAT(2,1)*BEAM(1)
     &      +OPTMAT(2,2)*BEAM(2)
     &      +OPTMAT(2,3)*BEAM(3)
     &      +OPTMAT(2,4)*BEAM(4)

          YBEAM=
     &      OPTMAT(3,1)*BEAM(1)
     &      +OPTMAT(3,2)*BEAM(2)
     &      +OPTMAT(3,3)*BEAM(3)
     &      +OPTMAT(3,4)*BEAM(4)

          TANTHEB=
     &      OPTMAT(4,1)*BEAM(1)
     &      +OPTMAT(4,2)*BEAM(2)
     &      +OPTMAT(4,3)*BEAM(3)
     &      +OPTMAT(4,4)*BEAM(4)

        ENDDO   !NELEM

        IF (ZBEAM*YBEAM.NE.0.0d0) THEN
          FOCUS=(EPSBEAM*EPSBEAM)/DABS(ZBEAM*YBEAM)
        ELSE
          WRITE(6,*)
     &      '*** ERROR IN PHASE_OMP: IMAGE IS IN FOCAL PLANE, CHECK INPUT ***'
          WRITE(LUNGFO,*)
     &      '*** ERROR IN PHASE_OMP: IMAGE IS IN FOCAL PLANE, CHECK INPUT ***'
          STOP
        ENDIF

        DO IOBS=1,NOBSV

          XOBS=OBSV(1,IOBS)
          YOBS=OBSV(2,IOBS)
          ZOBS=OBSV(3,IOBS)

          DO ISOUR=1,NSOURCE

            XSOUR=SOURCEN(1,1,ISOUR)
            IF (XOBS.LE.XSOUR) THEN
              WRITE(LUNGFO,*)'*** ERROR IN PHASE_OMP: Bad PINCEN(1)   ***'
              WRITE(LUNGFO,*)'CHECK INPUT FILE'
              WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
              WRITE(6,*)'*** ERROR IN PHASE_OMP: Bad PINCEN(1)  ***'
              WRITE(6,*)'CHECK INPUT FILE'
              WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
              STOP
            ENDIF

            YSOUR=SOURCEN(2,1,ISOUR)
            ZSOUR=SOURCEN(3,1,ISOUR)

            DX=(XSOUR-XOBS)
            DY=(YSOUR-YOBS)
            DZ=(ZSOUR-ZOBS)

            TANPHI=DZ/DX
            TANTHE=DY/DX

            PHI=ATAN2(DZ,DX)
            THETA=ATAN2(DY,DX)

            XPH=PHCENX
            DXPH=(XPH-XOBS)
            YPH=YOBS+TANTHE*DXPH
            ZPH=ZOBS+TANPHI*DXPH

            DR2PH=(XSOUR-XPH)**2+(YSOUR-YPH)**2+(ZSOUR-ZPH)**2
            DR2SOUR=DX**2+DY**2+DZ**2

            DO ifrq=1,NFREQ

              TGEO(1)=XPH
              TGEO(2)=YPH
              TGEO(3)=ZPH
              if(abs(tgeo(2)).lt.1.0d-15) tgeo(2)=0.0d0
              if(abs(tgeo(3)).lt.1.0d-15) tgeo(3)=0.0d0
              TGEO(4)=-TANTHE
              TGEO(5)=-TANPHI
              TGEO(6)=FREQ(ifrq)
              TGEO(7)=ifrq
              TGEO(8)=ISOUR
              TGEO(9)=XSOUR
              TGEO(10)=YSOUR
              TGEO(11)=ZSOUR
              if(abs(tgeo(10)).lt.1.0d-15) tgeo(10)=0.0d0
              if(abs(tgeo(11)).lt.1.0d-15) tgeo(11)=0.0d0
              ILIOBFR=ISOUR+NSOURCE*(IOBS-1+NOBSV*(ifrq-1))
              TGEO(12)=SPEC(ILIOBFR)
     &          *DR2SOUR/DR2PH
              TGEO(13)=XOBS
              TGEO(14)=YOBS
              TGEO(15)=ZOBS
              TGEO(16)=SPEC(ILIOBFR)
              CALL hfm(NIDGEO,TGEO)

              PHGEOSUM(ifrq)=PHGEOSUM(ifrq)+SPEC(ILIOBFR)*DA
            ENDDO !ifrq=1,NFREQ

            IF (
     &          TANTHE.LE.PHAPERYPP
     &          .AND.
     &          TANTHE.GE.PHAPERYPM
     &          .AND.
     &          TANPHI.LE.PHAPERZPP
     &          .AND.
     &          TANPHI.GE.PHAPERZPM
     &          .AND.
     &          YPH.LE.PHAPERYP
     &          .AND.
     &          YPH.GE.PHAPERYM
     &          .AND.
     &          ZPH.LE.PHAPERZP
     &          .AND.
     &          ZPH.GE.PHAPERZM
     &          ) THEN

              DO ifrq=1,NFREQ

                TGEO(1)=XPH
                TGEO(2)=YPH
                TGEO(3)=ZPH
                if(abs(tgeo(2)).lt.1.0d-15) tgeo(2)=0.0d0
                if(abs(tgeo(3)).lt.1.0d-15) tgeo(3)=0.0d0
                TGEO(4)=-TANTHE
                TGEO(5)=-TANPHI
                TGEO(6)=FREQ(ifrq)
                TGEO(7)=ifrq
                TGEO(8)=ISOUR
                TGEO(9)=XSOUR
                TGEO(10)=YSOUR
                TGEO(11)=ZSOUR
                if(abs(tgeo(10)).lt.1.0d-15) tgeo(10)=0.0d0
                if(abs(tgeo(11)).lt.1.0d-15) tgeo(11)=0.0d0
                ILIOBFR=ISOUR+NSOURCE*(IOBS-1+NOBSV*(ifrq-1))
                TGEO(12)=SPEC(ILIOBFR)
     &            *DR2SOUR/DR2PH
                TGEO(13)=XOBS
                TGEO(14)=YOBS
                TGEO(15)=ZOBS
                TGEO(16)=SPEC(ILIOBFR)
                CALL hfm(NIDGEO1,TGEO)

                PHGEOSEL(ifrq)=PHGEOSEL(ifrq)+SPEC(ILIOBFR)*DA
                ILIFR=ISOUR+NSOURCE*(ifrq-1)
                W=SPEC(ILIOBFR)
                PHMEANZ(ILIFR)=PHMEANZ(ILIFR)
     &            +ZPH*W
                PHSIGZ(ILIFR)=PHSIGZ(ILIFR)
     &            +ZPH*ZPH*W
                PHMEANY(ILIFR)=PHMEANY(ILIFR)
     &            +YPH*W
                PHSIGY(ILIFR)=PHSIGY(ILIFR)
     &            +YPH*YPH*W
                WSUM(ILIFR)=WSUM(ILIFR)+W

              ENDDO   !ifrq

C--- APPLY MATRICES OF BEAMLINE

              XBEAM=XPH
              YBEAM=YPH
              ZBEAM=ZPH
              TANTHEB=-TANTHE
              TANPHIB=-TANPHI

              DO J=1,4
                DO I=1,4
                  OPTMAT(I,J)=0.
                ENDDO
              ENDDO

              DO IELEM=1,NPHELEM

                BEAM(1)=ZBEAM
                BEAM(2)=TANPHIB
                BEAM(3)=YBEAM
                BEAM(4)=TANTHEB
                DO J=1,4
                  DO I=1,4
                    OPTMAT(I,J)=PHELEM(I,J,IELEM)
                  ENDDO
                ENDDO

                IF (
     &              BEAM(1).LT.PHELEM(5,1,IELEM)
     &              .OR.
     &              BEAM(1).GT.PHELEM(5,2,IELEM)
     &              .OR.
     &              BEAM(2).LT.PHELEM(5,3,IELEM)
     &              .OR.
     &              BEAM(2).GT.PHELEM(5,4,IELEM)
     &              ) THEN
                  GOTO 90   !OUT OF APERTURE, SKIP BEAM
                ELSE    !APERTURE CUT

                  ZBEAM=
     &              OPTMAT(1,1)*BEAM(1)
     &              +OPTMAT(1,2)*BEAM(2)
     &              +OPTMAT(1,3)*BEAM(3)
     &              +OPTMAT(1,4)*BEAM(4)

                  TANPHIB=
     &              OPTMAT(2,1)*BEAM(1)
     &              +OPTMAT(2,2)*BEAM(2)
     &              +OPTMAT(2,3)*BEAM(3)
     &              +OPTMAT(2,4)*BEAM(4)

                  YBEAM=
     &              OPTMAT(3,1)*BEAM(1)
     &              +OPTMAT(3,2)*BEAM(2)
     &              +OPTMAT(3,3)*BEAM(3)
     &              +OPTMAT(3,4)*BEAM(4)

                  TANTHEB=
     &              OPTMAT(4,1)*BEAM(1)
     &              +OPTMAT(4,2)*BEAM(2)
     &              +OPTMAT(4,3)*BEAM(3)
     &              +OPTMAT(4,4)*BEAM(4)

                  XBEAM=XBEAM+PHELEM(1,2,IELEM)

                ENDIF   !APERTURE CUT

              ENDDO   !NELEM

              DO ifrq=1,NFREQ

                TBEAM(1)=XBEAM
                TBEAM(2)=YBEAM
                TBEAM(3)=ZBEAM
                if(abs(tgeo(2)).lt.1.0d-15) tgeo(2)=0.0d0
                if(abs(tgeo(3)).lt.1.0d-15) tgeo(3)=0.0d0
                TBEAM(4)=TANTHEB
                TBEAM(5)=TANPHIB
                TBEAM(6)=FREQ(ifrq)
                TBEAM(7)=ifrq
                TBEAM(8)=ISOUR
                TBEAM(9)=XSOUR
                TBEAM(10)=YSOUR
                TBEAM(11)=ZSOUR
                if(abs(tgeo(10)).lt.1.0d-15) tgeo(10)=0.0d0
                if(abs(tgeo(11)).lt.1.0d-15) tgeo(11)=0.0d0
                ILIOBFR=ISOUR+NSOURCE*(IOBS-1+NOBSV*(ifrq-1))
                TBEAM(12)=SPEC(ILIOBFR)*DR2SOUR/DR2PH*FOCUS
                TBEAM(13)=XOBS
                TBEAM(14)=YOBS
                TBEAM(15)=ZOBS
                TBEAM(16)=SPEC(ILIOBFR)
                CALL hfm(NIDGEO2,TBEAM)
                PHBEAM(ifrq)=PHBEAM(ifrq)+SPEC(ILIOBFR)*DA
              ENDDO   !ifrq

90            CONTINUE

            ENDIF   !CUTS

          ENDDO  !ISOUR=1,NSOURCE

        ENDDO !IOBS=1,NOBSV

        CALL MHROUT(NIDGEO,ICYCLE,' ')
        CALL hdeletm(NIDGEO)
        CALL MHROUT(NIDGEO1,ICYCLE,' ')
        CALL hdeletm(NIDGEO1)
        CALL MHROUT(NIDGEO2,ICYCLE,' ')
        CALL hdeletm(NIDGEO2)

        DO ifrq=1,NFREQ
          DO ISOUR=1,NSOURCE
            ILIFR=ISOUR+NSOURCE*(ifrq-1)
            IF (WSUM(ILIFR).NE.0.0d0) THEN
              W=1.D0/WSUM(ILIFR)
              PHMEANZ(ILIFR)=PHMEANZ(ILIFR)*W
              PHMEANY(ILIFR)=PHMEANY(ILIFR)*W
              PHSIGZ(ILIFR)=DSQRT(DABS(
     &          PHSIGZ(ILIFR)*W
     &          -PHMEANZ(ILIFR)*PHMEANZ(ILIFR)))
              PHSIGY(ILIFR)=DSQRT(DABS(
     &          PHSIGY(ILIFR)*W
     &          -PHMEANY(ILIFR)*PHMEANY(ILIFR)))
            ENDIF   !(WSUM(ILIFR))
          ENDDO
        ENDDO

      ENDIF !ABS(IPHASE.GT.1)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      PHASE_OMP:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      PHCENX, PHCENY, PHCENZ:'
      WRITE(LUNGFO,*)'      ', PHCENX, PHCENY, PHCENZ
      WRITE(LUNGFO,*)'       PHWID, PHHIG:'
      WRITE(LUNGFO,*)'      ', PHWID, PHHIG
      WRITE(LUNGFO,*)'      NPHASEZ, NPHASEY:', NPHASEZ, NPHASEY
      WRITE(LUNGFO,*)'      MPHASEZ, MPHASEY:', MPHASEZ, MPHASEY
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      PHAPERZM,  PHAPERZP: ',PHAPERZM, PHAPERZP
      WRITE(LUNGFO,*)'      PHAPERYM,  PHAPERYP: ',PHAPERYM, PHAPERYP
      WRITE(LUNGFO,*)'      PHAPERZPM, PHAPERZPP:',PHAPERZPM, PHAPERZPP
      WRITE(LUNGFO,*)'      PHAPERYPM, PHAPERYPP:',PHAPERYPM, PHAPERYPP
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      IPHFOLD:         ',iphfold
      WRITE(LUNGFO,*)'      PHBETH, PHBETV:  ',PHBETH,PHBETV
      WRITE(LUNGFO,*)'      PHGSIGZ, PHGSIGY:',PHGSIGZ,PHGSIGY
      WRITE(LUNGFO,*)

      IF (ABS(IPHASE).GT.1) THEN

        WRITE(LUNGFO,*)'      NPHELEM: ',NPHELEM
        DO IELEM=1,NPHELEM
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'      ',(SNGL(PHELEM(1,I,IELEM)),I=1,4)
          WRITE(LUNGFO,*)'      ',(SNGL(PHELEM(2,I,IELEM)),I=1,4)
          WRITE(LUNGFO,*)'      ',(SNGL(PHELEM(3,I,IELEM)),I=1,4)
          WRITE(LUNGFO,*)'      ',(SNGL(PHELEM(4,I,IELEM)),I=1,4)
          WRITE(LUNGFO,*)'      ',(SNGL(PHELEM(5,I,IELEM)),I=1,4)
        ENDDO  !NPHELEM

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      focus factor:',SNGL(FOCUS)
        WRITE(LUNGFO,*)

        WRITE(LUNGFO,*)'      resulting matrix:'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      ',(SNGL(TOTMAT(1,I)),I=1,4)
        WRITE(LUNGFO,*)'      ',(SNGL(TOTMAT(2,I)),I=1,4)
        WRITE(LUNGFO,*)'      ',(SNGL(TOTMAT(3,I)),I=1,4)
        WRITE(LUNGFO,*)'      ',(SNGL(TOTMAT(4,I)),I=1,4)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '         photon energy, total intensity, selected intensity, ratio'
        WRITE(LUNGFO,*)
     &    '         (values refere to observation plane):'
        WRITE(LUNGFO,*)
        DO ifrq=1,NFREQ

          IF (PHGEOSUM(ifrq).NE.0.0d0) THEN
            SELGEO=PHGEOSEL(ifrq)/PHGEOSUM(ifrq)
          ELSE
            SELGEO=0.
          ENDIF
          WRITE(LUNGFO,*)'      ',SNGL(FREQ(ifrq))
     &      ,SNGL(PHGEOSUM(ifrq))
     &      ,SNGL(PHGEOSEL(ifrq))
     &      ,SELGEO
        ENDDO

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '         photon energy, selected intensity, beamline intensity, ratio'
        WRITE(LUNGFO,*)
     &    '         (values refere to image plane):'
        WRITE(LUNGFO,*)
        DO ifrq=1,NFREQ

          IF (PHGEOSEL(ifrq).NE.0.0d0) THEN
            SELGEO=PHBEAM(ifrq)/PHGEOSEL(ifrq)
          ELSE
            SELGEO=0.
          ENDIF
          WRITE(LUNGFO,*)'      ',SNGL(FREQ(ifrq))
     &      ,SNGL(PHGEOSEL(ifrq))
     &      ,SNGL(PHBEAM(ifrq))
     &      ,SELGEO
        ENDDO
        WRITE(LUNGFO,*)

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      Effective source size in phase plane:'
        WRITE(LUNGFO,*)'      ifrq     ISOURCE     Zrms     Yrms'
        WRITE(LUNGFO,*)
        DO ifrq=1,NFREQ
          DO ISOUR=1,NSOURCE
            ILIFR=ISOUR+NSOURCE*(ifrq-1)
            WRITE(LUNGFO,*)'      ',ifrq,ISOUR
     &        ,SNGL(PHSIGZ(ILIFR))
     &        ,SNGL(PHSIGY(ILIFR))
            TSIZ(1)=ifrq
            TSIZ(2)=ISOUR
            TSIZ(3)=PHSIGZ(ILIFR)
            TSIZ(4)=PHSIGY(ILIFR)
            CALL hfm(NIDGEO+3,TSIZ)
          ENDDO
        ENDDO
        CALL MHROUT(NIDGEO+3,ICYCLE,' ')
        CALL hdeletm(NIDGEO+3)

      ENDIF !ABS(IPHASE.GT.1)

      WRITE(LUNGFO,*)

      IF (IHSEL.NE.0) THEN

        IF (ifreq2P.EQ.2) THEN

          FLOW=DLOG10(FREQ(1)/1.5)
          FHIG=DLOG10(FREQ(NFREQ)*1.5)

          call hbook1m(IDSEL,'INTEGRATED PHASESPACE (GEOM. OPTIC)',
     &      NFREQ,FLOW,FHIG,0.)
          call hbook1m(IDSEL+1,
     &      'INTEGRATED ACCEPTED PHASESPACE (GEOM. OPTIC)',
     &      NFREQ,FLOW,FHIG,0.)
          call hbook1m(IDSEL+2,
     &      'RATIO OF INTEGRATED AND ACCEPTED PHASESPACE',
     &      NFREQ,FLOW,FHIG,0.)
          call hbook1m(IDSEL+10,
     &      'INTEGRATED ACCEPTED PHASESPACE AT END OF BEAMLINE',
     &      NFREQ,FLOW,FHIG,0.)
          call hbook1m(IDSEL+11,
     &      'RATIO OF INTEGRATED AND PHASESPACE AFTER BEAMLINE',
     &      NFREQ,FLOW,FHIG,0.)

          DO ifrq=1,NFREQ
            IF (PHGEOSUM(ifrq).GT.0.0d0)
     &        CALL hfillm
     &        (IDSEL,SNGL(FREQ(ifrq)),0.,DLOG10(PHGEOSUM(ifrq)))
            IF (PHGEOSEL(ifrq).GT.0.0d0)
     &        CALL hfillm
     &        (IDSEL+1,SNGL(FREQ(ifrq)),0.,DLOG10(PHGEOSEL(ifrq)))
            IF (PHBEAM(ifrq).GT.0.0d0)
     &        CALL hfillm
     &        (IDSEL+10,SNGL(FREQ(ifrq)),0.,DLOG10(PHBEAM(ifrq)))
          ENDDO

        ELSE

          DF=FREQ(2)-FREQ(1)
          IF (DF.EQ.0.) DF=1.

          FLOW=FREQ(1)-DF/2.
          FHIG=FREQ(NFREQ)+DF/2.

          call hbook1m(IDSEL,'INTEGRATED PHASESPACE (GEOM. OPTIC)',
     &      NFREQ,FLOW,FHIG,0.)
          call hbook1m(IDSEL+1,
     &      'INTEGRATED ACCEPTED PHASESPACE (GEOM. OPTIC)',
     &      NFREQ,FLOW,FHIG,0.)
          call hbook1m(IDSEL+2,
     &      'RATIO OF INTEGRATED AND ACCEPTED PHASESPACE',
     &      NFREQ,FLOW,FHIG,0.)
          call hbook1m(IDSEL+10,
     &      'INTEGRATED ACCEPTED PHASESPACE AT END OF BEAMLINE',
     &      NFREQ,FLOW,FHIG,0.)
          call hbook1m(IDSEL+11,
     &      'RATIO OF INTEGRATED AND PHASESPACE AFTER BEAMLINE',
     &      NFREQ,FLOW,FHIG,0.)

          DO ifrq=1,NFREQ
            CALL hfillm
     &        (IDSEL,SNGL(FREQ(ifrq)),0.,PHGEOSUM(ifrq))
            CALL hfillm
     &        (IDSEL+1,SNGL(FREQ(ifrq)),0.,PHGEOSEL(ifrq))
            CALL hfillm
     &        (IDSEL+10,SNGL(FREQ(ifrq)),0.,PHBEAM(ifrq))
          ENDDO

        ENDIF  !ifreq2P

        CALL hoperam(IDSEL+1,'/',IDSEL,IDSEL+2,1.,1.)
        CALL hoperam(IDSEL+10,'/',IDSEL,IDSEL+11,1.,1.)

        CALL MHROUT(IDSEL,ICYCLE,' ')
        CALL MHROUT(IDSEL+1,ICYCLE,' ')
        CALL MHROUT(IDSEL+2,ICYCLE,' ')
        CALL MHROUT(IDSEL+10,ICYCLE,' ')
        CALL MHROUT(IDSEL+11,ICYCLE,' ')
        CALL hdeletm(IDSEL)
        CALL hdeletm(IDSEL+1)
        CALL hdeletm(IDSEL+2)
        CALL hdeletm(IDSEL+10)
        CALL hdeletm(IDSEL+11)

      ENDIF !IHSEL.NE.0

      if (mhbookp.eq.0.and.iroottrees.ge.0) then
        CALL hrendm('PHASE')
        CLOSE(LUNPH)
        CALL hcdirm(OLDDIR,' ')
      endif

      DEALLOCATE(WSUM)
      DEALLOCATE(PHMEANZ)
      DEALLOCATE(PHMEANY)
      DEALLOCATE(PHSIGZ)
      DEALLOCATE(PHSIGY)
      DEALLOCATE(PHSHIFT)
      DEALLOCATE(AMPLI)
      DEALLOCATE(phspec3)
      DEALLOCATE(phspec3f,phspec3fy)
      deALLOCATE(specwz)
      deALLOCATE(specfwz)
      deALLOCATE(yphw)
      deALLOCATE(specwy)
      deALLOCATE(specfwy)
      DEALLOCATE(EXPOM)
      DEALLOCATE(DEXPOM)
      DEALLOCATE(phws1)
      DEALLOCATE(phws2)
      DEALLOCATE(phws3)
      DEALLOCATE(phws4)
      DEALLOCATE(phcoef)
      DEALLOCATE(freq_omp)

      mthreads=mthreadso

      RETURN
      END
