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
      SUBROUTINE PHASE
*KEEP,gplhint.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
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

* ROUTINE TO PROPAGATE COMPLEXE AMPLITUDE FROM PINHOLE BACK TO
* LOCATION OF EFFECTIVE SOURCE AT (PHCENX,PHCENY,PHCENZ)


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
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

      CHARACTER(8) OLDDIR

      INTEGER ICYCLE,NTUP_P,IOBS,IPHZ,IPHY,IFREQ,IEPS,I,NGEO_P,ISOUR
      INTEGER NIDGEO1,ISTAT,NIDGEO2,NBEAM_P,J,IELEM,NSIZE_P
      INTEGER IOBSY,IOBSZ,K,ix

      PARAMETER(NTUP_P=15,NGEO_P=16,NBEAM_P=16,NSIZE_P=4)

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
     &        ,YAI2(NDOBSVZP)
     &  ,YAI3(NDOBSVZP)
     &  ,RESULT(3,2)

      DOUBLE PRECISION XAY(NDOBSVYP)
     &  ,YAR1Y(NDOBSVYP)
     &  ,YAR2Y(NDOBSVYP)
     &  ,YAR3Y(NDOBSVYP)
     &  ,YAI1Y(NDOBSVYP)
     &  ,YAI2Y(NDOBSVYP)
     &  ,YAI3Y(NDOBSVYP)
     &  ,RESULTY(3,2)

      DOUBLE PRECISION WS1(NDOBSVZP+NDOBSVYP),WS2(NDOBSVZP+NDOBSVYP)
      DOUBLE PRECISION WS3(NDOBSVZP+NDOBSVYP),WS4(NDOBSVZP+NDOBSVYP)
      DOUBLE PRECISION COEFF(NDOBSVZP+NDOBSVYP)

      double precision, dimension (:), allocatable :: zphw,yphw,
     &  specwz,specwy,specfwz,specfwy,phws1,phws2,phcoef,phws3,phws4
      double precision phgsigz,phgsigy,wlen,sigrp,sigr

      integer is0

      character(8) chphase
      integer lenchphase

      DATA CHTAGS
     &  /'x','y','z','e','ie','iy','iz',
     &  're_x','im_x','re_y','im_y','re_z','im_z','spec','specf'/
      data chgeo
     &  /'x','y','z','yp','zp','e','ie','is','xs','ys','zs','spec'
     &  ,'xo','yo','zo','speco'/
      data chbeam
     &  /'x','y','z','yp','zp','e','ie','is','xs','ys','zs','spec'
     &  ,'xo','yo','zo','speco'/
      data chsize /'ie','is','zrms','yrms'/

      DATA EPSBEAM/0.001D0/


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
      sigrp=sqrt(wlen/(sourceeo(1,1,1)-sourceao(1,1,1)))
      sigr=wlen/twopi1/sigrp
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

      DA=PINW/(MOBSVZ-1)*PINH/(MOBSVY-1)

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
      ALLOCATE(AMPLI(3,mphasez,mphasey,NFREQ))
      ALLOCATE(phspec(mphasez,mphasey))
      ALLOCATE(phspecf(mphasez,mphasey))
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

      ampli=(0.0d0,0.0d0)
      specwz=0.0d0
      specwy=0.0d0
      phspec=0.0d0
      phspecf=0.0d0

      if (mhbookp.eq.0 .and. iroottrees.ge.0) then
        CALL hcdirm(OLDDIR,'R')
        CALL hropenm(LUNPH,'PHASE',FILEPH,'N',1024,ISTAT)
        CALL hcdirm(chphase(1:lenchphase),' ')
        IF (ISTAT.NE.0) THEN
          WRITE(6,*)'*** ERROR IN hropenm (SR PHASE) ***'
          WRITE(LUNGFO,*)'*** ERROR IN hropenm (SR PHASE) ***'
          STOP
        ENDIF
        CALL MHROUT(IDCODE,ICYCLE,' ')
      endif

      DO IFREQ=1,NFREQ
        PHGEOSUM(IFREQ)=0.D0
        PHGEOSEL(IFREQ)=0.D0
        PHBEAM(IFREQ)=0.D0
        DO ISOUR=1,NSOURCE
          ILIFR=ISOUR+NSOURCE*(IFREQ-1)
          WSUM(ILIFR)=0.D0
          PHMEANZ(ILIFR)=0.D0
          PHMEANY(ILIFR)=0.D0
          PHSIGZ(ILIFR)=0.D0
          PHSIGY(ILIFR)=0.D0
        ENDDO
      ENDDO

      CALL hbookm(NIDPHASE,'PHASE',NTUP_P,chphase(1:lenchphase),
     &  mphasez*mphasey*nfreq,CHTAGS)
      XPH=PHCENX
      XOBS=PINCEN(1)
      DX=XOBS-XPH
      DX2=DX*DX

      IF (DX.EQ.0.0D0) THEN
        WRITE(LUNGFO,*)'*** ERROR IN PHASE: PHCENX=PINCEN(1)  ***'
        WRITE(LUNGFO,*)'CHECK INPUT FILE'
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(6,*)'*** ERROR IN PHASE: PHCENX=PINCEN(1)  ***'
        WRITE(6,*)'CHECK INPUT FILE'
        WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
        STOP
      ENDIF

      PHLOWZ=PHCENZ-PHWID/2.D0
      PHLOWY=PHCENY-PHHIG/2.D0
      YPH=PHLOWY-DMASHY

      OMC=FREQ(1)/(HBAREV1*CLIGHT1)
      IF (IFREQ2P.GT.2) THEN
        DOMC=(FREQ(2)-FREQ(1))/(HBAREV1*CLIGHT1)
      ELSE
        DOMC=OMC
      ENDIF !(IFREQ2P.GT.2)

      DO iphy=(mphasey-nphasey)/2+1,(mphasey-nphasey)/2+NPHASEY

        YPH=YPH+DMASHY
        ZPH=PHLOWZ-DMASHZ

        DO IPHZ=(mphasez-nphasez)/2+1,nphasez+(mphasez-nphasez)/2

          ZPH=ZPH+DMASHZ

          DO IOBS=1,NOBSV

            XOBS=OBSV(1,IOBS)
            YOBS=OBSV(2,IOBS)
            ZOBS=OBSV(3,IOBS)

            DY=YOBS-YPH
            DZ=ZOBS-ZPH
            DZY2=DZ*DZ+DY*DY

C     TO MAKE SURE THAT TAYLOR-EXPANSION IS VALID

            IF (DZY2.GT.0.01D0*DX2) THEN
              WRITE(LUNGFO,*)'*** ERROR IN PHASE: DZY2.GT.0.01D0*DX2  ***'
              WRITE(LUNGFO,*)'CHECK INPUT FILE AND INCREASE PINCEN(1)'
              WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
              WRITE(6,*)'*** ERROR IN PHASE: PHCENX=PINCEN(1)  ***'
              WRITE(6,*)'CHECK INPUT FILE AND INCREASE PINCEN(1)'
              WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
              STOP
            ENDIF

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
     &        -0.0390625D0*eps(4)+
     &        0.0625D0*eps(3)-0.125D0*eps(2)+0.5D0*eps(1)

            DR=DABS(DX*(ANS+1.D0))
            DRRED=-DABS(DX*ANS)

            IF (DR.NE.0.D0) THEN
              EXPOM(IOBS)=CDEXP(DCMPLX(0.D0,DRRED*OMC))/DR
            ELSE
              EXPOM(IOBS)=1.D0
            ENDIF
            DEXPOM(IOBS)=CDEXP(DCMPLX(0.D0,DRRED*DOMC))

          ENDDO   !NOBS

          DO IFREQ=1,NFREQ

            RLAMBDA1=FREQ(IFREQ)/WTOE1*1.D9   !1/lambda[m]=1/(wtoe1/freq*1.e-9)

            IF (IPHASE.GT.0) THEN

              DO IOBS=1,NOBSV

                IOBFR=IOBS+NOBSV*(IFREQ-1)

                IF (IFREQ.EQ.1) THEN

                  PHSHIFT(IOBS)=EXPOM(IOBFR)
                ELSE
                  PHSHIFT(IOBS)=PHSHIFT(IOBS)*DEXPOM(IOBS)
                ENDIF   !(IFREQ.EQ.1)


                IF (DX.GE.0) THEN

                  ampli(1,iphz,iphy,IFREQ)=ampli(1,iphz,iphy,IFREQ)+
     &              DCMPLX(REAIMA(1,1,IOBFR),REAIMA(1,2,IOBFR))
     &              *PHSHIFT(IOBS)
                  ampli(2,iphz,iphy,IFREQ)=ampli(2,iphz,iphy,IFREQ)+
     &              DCMPLX(REAIMA(2,1,IOBFR),REAIMA(2,2,IOBFR))
     &              *PHSHIFT(IOBS)
                  ampli(3,iphz,iphy,IFREQ)=ampli(3,iphz,iphy,IFREQ)+
     &              DCMPLX(REAIMA(3,1,IOBFR),REAIMA(3,2,IOBFR))
     &              *PHSHIFT(IOBS)

                ELSE

                  ampli(1,iphz,iphy,IFREQ)=ampli(1,iphz,iphy,IFREQ)+
     &              DCMPLX(REAIMA(1,1,IOBFR),-REAIMA(1,2,IOBFR))
     &              *PHSHIFT(IOBS)
                  ampli(2,iphz,iphy,IFREQ)=ampli(2,iphz,iphy,IFREQ)+
     &              DCMPLX(REAIMA(2,1,IOBFR),-REAIMA(2,2,IOBFR))
     &              *PHSHIFT(IOBS)
                  ampli(3,iphz,iphy,IFREQ)=ampli(3,iphz,iphy,IFREQ)+
     &              DCMPLX(REAIMA(3,1,IOBFR),-REAIMA(3,2,IOBFR))
     &              *PHSHIFT(IOBS)

                ENDIF !(DX.GE.0)


              ENDDO  !NOBSV

              ampli(1,iphz,iphy,IFREQ)=ampli(1,iphz,iphy,IFREQ)*DA*RLAMBDA1
              ampli(2,iphz,iphy,IFREQ)=ampli(2,iphz,iphy,IFREQ)*DA*RLAMBDA1
              ampli(3,iphz,iphy,IFREQ)=ampli(3,iphz,iphy,IFREQ)*DA*RLAMBDA1

            ELSE  !IPHASE.GT.0

              IOBS=0
              DO IOBSY=1,NOBSVY
                DO IOBSZ=1,NOBSVZ

                  IOBS=IOBS+1

                  IF (IFREQ.EQ.1) THEN
                    PHSHIFT(IOBS)=EXPOM(IOBS+NOBSV*(IFREQ-1))
                  ELSE
                    PHSHIFT(IOBS)=PHSHIFT(IOBS)*DEXPOM(IOBS)
                  ENDIF !(IFREQ.EQ.1)

                  IF (DX.GE.0) THEN

                    IOBFR=IOBS+NOBSV*(IFREQ-1)
                    ampli(1,iphz,iphy,IFREQ)=
     &                DCMPLX(REAIMA(1,1,IOBFR),REAIMA(1,2,IOBFR))
     &                *PHSHIFT(IOBS)
                    ampli(2,iphz,iphy,IFREQ)=
     &                DCMPLX(REAIMA(2,1,IOBFR),REAIMA(2,2,IOBFR))
     &                *PHSHIFT(IOBS)
                    ampli(3,iphz,iphy,IFREQ)=
     &                DCMPLX(REAIMA(3,1,IOBFR),REAIMA(3,2,IOBFR))
     &                *PHSHIFT(IOBS)

                  ELSE

                    ampli(1,iphz,iphy,IFREQ)=
     &                DCMPLX(REAIMA(1,1,IOBFR),-REAIMA(1,2,IOBFR))
     &                *PHSHIFT(IOBS)
                    ampli(2,iphz,iphy,IFREQ)=
     &                DCMPLX(REAIMA(2,1,IOBFR),-REAIMA(2,2,IOBFR))
     &                *PHSHIFT(IOBS)
                    ampli(3,iphz,iphy,IFREQ)=
     &                DCMPLX(REAIMA(3,1,IOBFR),-REAIMA(3,2,IOBFR))
     &                *PHSHIFT(IOBS)

                  ENDIF !(DX.GE.0)

                  XA(IOBSZ)=OBSVZ(IOBSZ)

                  YAR1(IOBSZ)=DREAL(ampli(1,iphz,iphy,IFREQ))
                  YAR2(IOBSZ)=DREAL(ampli(2,iphz,iphy,IFREQ))
                  YAR3(IOBSZ)=DREAL(ampli(3,iphz,iphy,IFREQ))

                  YAI1(IOBSZ)=DIMAG(ampli(1,iphz,iphy,IFREQ))
                  YAI2(IOBSZ)=DIMAG(ampli(2,iphz,iphy,IFREQ))
                  YAI3(IOBSZ)=DIMAG(ampli(3,iphz,iphy,IFREQ))

                ENDDO   !NOBSVZ

                CALL UTIL_SPLINE_INTEGRAL
     &            (XA,YAR1,NOBSVZ,RESULT(1,1),COEFF,WS1,WS2,WS3,WS4)
                CALL UTIL_SPLINE_INTEGRAL
     &            (XA,YAI1,NOBSVZ,RESULT(1,2),COEFF,WS1,WS2,WS3,WS4)

                CALL UTIL_SPLINE_INTEGRAL
     &            (XA,YAR2,NOBSVZ,RESULT(2,1),COEFF,WS1,WS2,WS3,WS4)
                CALL UTIL_SPLINE_INTEGRAL
     &            (XA,YAI2,NOBSVZ,RESULT(2,2),COEFF,WS1,WS2,WS3,WS4)

                CALL UTIL_SPLINE_INTEGRAL
     &            (XA,YAR3,NOBSVZ,RESULT(3,1),COEFF,WS1,WS2,WS3,WS4)
                CALL UTIL_SPLINE_INTEGRAL
     &            (XA,YAI3,NOBSVZ,RESULT(3,2),COEFF,WS1,WS2,WS3,WS4)

                XAY(IOBSY)=OBSVY(IOBSY)
                YAR1Y(IOBSY)=RESULT(1,1)
                YAI1Y(IOBSY)=RESULT(1,2)
                YAR2Y(IOBSY)=RESULT(2,1)
                YAI2Y(IOBSY)=RESULT(2,2)
                YAR3Y(IOBSY)=RESULT(3,1)
                YAI3Y(IOBSY)=RESULT(3,2)

              ENDDO !NOBSVY

              CALL UTIL_SPLINE_INTEGRAL
     &          (XAY,YAR1Y,NOBSVY,RESULTY(1,1),COEFF,WS1,WS2,WS3,WS4)
              CALL UTIL_SPLINE_INTEGRAL
     &          (XAY,YAI1Y,NOBSVY,RESULTY(1,2),COEFF,WS1,WS2,WS3,WS4)

              CALL UTIL_SPLINE_INTEGRAL
     &          (XAY,YAR2Y,NOBSVY,RESULTY(2,1),COEFF,WS1,WS2,WS3,WS4)
              CALL UTIL_SPLINE_INTEGRAL
     &          (XAY,YAI2Y,NOBSVY,RESULTY(2,2),COEFF,WS1,WS2,WS3,WS4)

              CALL UTIL_SPLINE_INTEGRAL
     &          (XAY,YAR3Y,NOBSVY,RESULTY(3,1),COEFF,WS1,WS2,WS3,WS4)
              CALL UTIL_SPLINE_INTEGRAL
     &          (XAY,YAI3Y,NOBSVY,RESULTY(3,2),COEFF,WS1,WS2,WS3,WS4)


              ampli(1,iphz,iphy,IFREQ)=DCMPLX(RESULTY(1,1),RESULTY(1,2))*RLAMBDA1
              ampli(2,iphz,iphy,IFREQ)=DCMPLX(RESULTY(2,1),RESULTY(2,2))*RLAMBDA1
              ampli(3,iphz,iphy,IFREQ)=DCMPLX(RESULTY(3,1),RESULTY(3,2))*RLAMBDA1

            ENDIF !IPHASE.GT.0

          ENDDO   !NFREQ

        ENDDO  !NPHASEZ
      ENDDO !NPHAZEY

      do iphz=1,mphasez
        zphw(iphz)=-(mphasez-1)*dmashz/2.0+(iphz-1)*dmashz
      enddo
      do iphy=1,mphasey
        yphw(iphy)=-(mphasey-1)*dmashy/2.0+(iphy-1)*dmashy
      enddo

      smax=-1.0d30
      do ifreq=1,nfreq

        YPH=PHLOWY-DMASHY
        DO iphy=(mphasey-nphasey)/2+1,(mphasey-nphasey)/2+NPHASEY

          YPH=YPH+DMASHY
          ZPH=PHLOWZ-DMASHZ

          DO IPHZ=(mphasez-nphasez)/2+1,(mphasez-nphasez)/2+NPHASEZ
            ZPH=ZPH+DMASHZ
            PHSPEC(iphz,iphy)=
     &        DREAL(ampli(1,iphz,iphy,IFREQ))*DREAL(ampli(1,iphz,iphy,IFREQ))+
     &        DIMAG(ampli(1,iphz,iphy,IFREQ))*DIMAG(ampli(1,iphz,iphy,IFREQ))+
     &        DREAL(ampli(2,iphz,iphy,IFREQ))*DREAL(ampli(2,iphz,iphy,IFREQ))+
     &        DIMAG(ampli(2,iphz,iphy,IFREQ))*DIMAG(ampli(2,iphz,iphy,IFREQ))+
     &        DREAL(ampli(3,iphz,iphy,IFREQ))*DREAL(ampli(3,iphz,iphy,IFREQ))+
     &        DIMAG(ampli(3,iphz,iphy,IFREQ))*DIMAG(ampli(3,iphz,iphy,IFREQ))

            if (phspec(iphz,iphy).gt.smax) smax=phspec(iphz,iphy)

          ENDDO   !NPHASEZ

        ENDDO  !NPHAZEY

        if (iphfold.ne.0) then

          do iphy=(mphasey-nphasey)/2+1,(mphasey-nphasey)/2+nphasey
c          do iphy=1,mphasey
            specwz=phspec(1:mphasez,iphy)
            if (dgsigz(1).gt.0.0d0.and.phgsigz.gt.0.0d0) then
              if (iphfold.lt.0) then
                call util_fold_function_gauss_lin(
     &            mphasez,zphw,specwz,phgsigz,dgsigz(1),specfwz,phws1,phws2)
              else
                call util_fold_function_gauss(
     &            mphasez,zphw,specwz,phgsigz,dgsigz(1),specfwz,
     &            phcoef,phws1,phws2,phws3,phws4)
              endif
            else
              specfwz=specwz
            endif
            phspecf(1:mphasez,iphy)=specfwz
          enddo !iphy

c          do iphz=(mphasez-nphasez)/2+1,(mphasez-nphasez)/2+nphasez
          do iphz=1,mphasez
            specwy=phspecf(iphz,1:mphasey)
            if (dgsigy(1).gt.0.0d0.and.phgsigy.gt.0.0d0) then
              if (iphfold.lt.0) then
                call util_fold_function_gauss_lin(
     &            mphasey,yphw,specwy,phgsigy,dgsigy(1),specfwy,phws1,phws2)
              else
                call util_fold_function_gauss(
     &            mphasey,yphw,specwy,phgsigy,dgsigy(1),specfwy,phcoef,
     &            phws1,phws2,phws3,phws4)
              endif
            else
              specfwy=specwy
            endif
            phspecf(iphz,1:mphasey)=specfwy
          enddo !iphz

c          phspecf=phspecf*specnor

        endif !iphfold

c        print*,'PHASE, SPECNOR: ',SPECNOR
c        phspec=phspec*specnor
c        smax=smax*specnor

c        YPH=PHLOWY-DMASHY
c        DO iphy=(mphasey-nphasey)/2+1,(mphasey-nphasey)/2+NPHASEY

        sfmax=-1.0d30

        DO iphy=1,mphasey

c          YPH=YPH+DMASHY
c          ZPH=PHLOWZ-DMASHZ

c          DO IPHZ=(mphasez-nphasez)/2+1,(mphasez-nphasez)/2+NPHASEZ
          DO IPHZ=1,mphasez
            zPH=zPH+DMASHz
            TUP(1)=XPH
c            TUP(2)=YPH
c            TUP(3)=ZPH
            TUP(2)=YPHw(iphy)
            TUP(3)=ZPHw(iphz)
            TUP(4)=FREQ(IFREQ)
            TUP(5)=IFREQ
            TUP(6)=IPHY
            TUP(7)=IPHZ
            TUP(8)=DREAL(ampli(1,iphz,iphy,IFREQ))
            TUP(9)=DIMAG(ampli(1,iphz,iphy,IFREQ))
            TUP(10)=DREAL(ampli(2,iphz,iphy,IFREQ))
            TUP(11)=DIMAG(ampli(2,iphz,iphy,IFREQ))
            TUP(12)=DREAL(ampli(3,iphz,iphy,IFREQ))
            TUP(13)=DIMAG(ampli(3,iphz,iphy,IFREQ))
            TUP(14)=PHSPEC(iphz,iphy)
            TUP(15)=PHSPECf(iphz,iphy)
            if (phspecf(iphz,iphy).gt.sfmax) sfmax=phspecf(iphz,iphy)
            CALL hfm(NIDPHASE,TUP)
          ENDDO   !mPHASEZ
        ENDDO  !MPHASEY

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
          WRITE(6,*)'*** ERROR IN hropenm (SR PHASE) ***'
          WRITE(LUNGFO,*)'*** ERROR IN hropenm (SR PHASE) ***'
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

        IF (ZBEAM*YBEAM.NE.0.D0) THEN
          FOCUS=(EPSBEAM*EPSBEAM)/DABS(ZBEAM*YBEAM)
        ELSE
          WRITE(6,*)
     &      '*** ERROR IN PHASE: IMAGE IS IN FOCAL PLANE, CHECK INPUT ***'
          WRITE(LUNGFO,*)
     &      '*** ERROR IN PHASE: IMAGE IS IN FOCAL PLANE, CHECK INPUT ***'
          STOP
        ENDIF

        DO IOBS=1,NOBSV

          XOBS=OBSV(1,IOBS)
          YOBS=OBSV(2,IOBS)
          ZOBS=OBSV(3,IOBS)

          DO ISOUR=1,NSOURCE

            XSOUR=SOURCEN(1,1,ISOUR)
            IF (XOBS.LE.XSOUR) THEN
              WRITE(LUNGFO,*)'*** ERROR IN PHASE: Bad PINCEN(1)   ***'
              WRITE(LUNGFO,*)'CHECK INPUT FILE'
              WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
              WRITE(6,*)'*** ERROR IN PHASE: Bad PINCEN(1)  ***'
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

            DO IFREQ=1,NFREQ

              TGEO(1)=XPH
              TGEO(2)=YPH
              TGEO(3)=ZPH
              TGEO(4)=-TANTHE
              TGEO(5)=-TANPHI
              TGEO(6)=FREQ(IFREQ)
              TGEO(7)=IFREQ
              TGEO(8)=ISOUR
              TGEO(9)=XSOUR
              TGEO(10)=YSOUR
              TGEO(11)=ZSOUR
              ILIOBFR=ISOUR+NSOURCE*(IOBS-1+NOBSV*(IFREQ-1))
              TGEO(12)=SPEC(ILIOBFR)
     &          *DR2SOUR/DR2PH
              TGEO(13)=XOBS
              TGEO(14)=YOBS
              TGEO(15)=ZOBS
              TGEO(16)=SPEC(ILIOBFR)
              CALL hfm(NIDGEO,TGEO)

              PHGEOSUM(IFREQ)=PHGEOSUM(IFREQ)+SPEC(ILIOBFR)*DA
            ENDDO !IFREQ=1,NFREQ

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

              DO IFREQ=1,NFREQ

                TGEO(1)=XPH
                TGEO(2)=YPH
                TGEO(3)=ZPH
                TGEO(4)=-TANTHE
                TGEO(5)=-TANPHI
                TGEO(6)=FREQ(IFREQ)
                TGEO(7)=IFREQ
                TGEO(8)=ISOUR
                TGEO(9)=XSOUR
                TGEO(10)=YSOUR
                TGEO(11)=ZSOUR
                ILIOBFR=ISOUR+NSOURCE*(IOBS-1+NOBSV*(IFREQ-1))
                TGEO(12)=SPEC(ILIOBFR)
     &            *DR2SOUR/DR2PH
                TGEO(13)=XOBS
                TGEO(14)=YOBS
                TGEO(15)=ZOBS
                TGEO(16)=SPEC(ILIOBFR)
                CALL hfm(NIDGEO1,TGEO)

                PHGEOSEL(IFREQ)=PHGEOSEL(IFREQ)+SPEC(ILIOBFR)*DA
                ILIFR=ISOUR+NSOURCE*(IFREQ-1)
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

              ENDDO   !IFREQ


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

              DO IFREQ=1,NFREQ

                TBEAM(1)=XBEAM
                TBEAM(2)=YBEAM
                TBEAM(3)=ZBEAM
                TBEAM(4)=TANTHEB
                TBEAM(5)=TANPHIB
                TBEAM(6)=FREQ(IFREQ)
                TBEAM(7)=IFREQ
                TBEAM(8)=ISOUR
                TBEAM(9)=XSOUR
                TBEAM(10)=YSOUR
                TBEAM(11)=ZSOUR
                ILIOBFR=ISOUR+NSOURCE*(IOBS-1+NOBSV*(IFREQ-1))
                TBEAM(12)=SPEC(ILIOBFR)*DR2SOUR/DR2PH*FOCUS
                TBEAM(13)=XOBS
                TBEAM(14)=YOBS
                TBEAM(15)=ZOBS
                TBEAM(16)=SPEC(ILIOBFR)
                CALL hfm(NIDGEO2,TBEAM)
                PHBEAM(IFREQ)=PHBEAM(IFREQ)+SPEC(ILIOBFR)*DA
              ENDDO   !IFREQ

90          CONTINUE

          ENDIF   !CUTS

        ENDDO  !ISOUR=1,NSOURCE

      ENDDO !IOBS=1,NOBSV

      CALL MHROUT(NIDGEO,ICYCLE,' ')
      CALL hdeletm(NIDGEO)
      CALL MHROUT(NIDGEO1,ICYCLE,' ')
      CALL hdeletm(NIDGEO1)
      CALL MHROUT(NIDGEO2,ICYCLE,' ')
      CALL hdeletm(NIDGEO2)

      DO IFREQ=1,NFREQ
        DO ISOUR=1,NSOURCE
            ILIFR=ISOUR+NSOURCE*(IFREQ-1)
            IF (WSUM(ILIFR).NE.0.D0) THEN
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
      WRITE(LUNGFO,*)'      SR PHASE:'
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
        DO IFREQ=1,NFREQ

          IF (PHGEOSUM(IFREQ).NE.0.D0) THEN
            SELGEO=PHGEOSEL(IFREQ)/PHGEOSUM(IFREQ)
          ELSE
            SELGEO=0.
          ENDIF
          WRITE(LUNGFO,*)'      ',SNGL(FREQ(IFREQ))
     &      ,SNGL(PHGEOSUM(IFREQ))
     &      ,SNGL(PHGEOSEL(IFREQ))
     &      ,SELGEO
        ENDDO


        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '         photon energy, selected intensity, beamline intensity, ratio'
        WRITE(LUNGFO,*)
     &    '         (values refere to image plane):'
        WRITE(LUNGFO,*)
        DO IFREQ=1,NFREQ

          IF (PHGEOSEL(IFREQ).NE.0.D0) THEN
            SELGEO=PHBEAM(IFREQ)/PHGEOSEL(IFREQ)
          ELSE
            SELGEO=0.
          ENDIF
          WRITE(LUNGFO,*)'      ',SNGL(FREQ(IFREQ))
     &      ,SNGL(PHGEOSEL(IFREQ))
     &      ,SNGL(PHBEAM(IFREQ))
     &      ,SELGEO
        ENDDO
        WRITE(LUNGFO,*)

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      Effective source size in phase plane:'
        WRITE(LUNGFO,*)'      IFREQ     ISOURCE     Zrms     Yrms'
        WRITE(LUNGFO,*)
        DO IFREQ=1,NFREQ
          DO ISOUR=1,NSOURCE
            ILIFR=ISOUR+NSOURCE*(IFREQ-1)
            WRITE(LUNGFO,*)'      ',IFREQ,ISOUR
     &        ,SNGL(PHSIGZ(ILIFR))
     &        ,SNGL(PHSIGY(ILIFR))
            TSIZ(1)=IFREQ
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

        IF (IFREQ2P.EQ.2) THEN

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

          DO IFREQ=1,NFREQ
            IF (PHGEOSUM(IFREQ).GT.0.D0)
     &        CALL hfillm
     &        (IDSEL,SNGL(FREQ(IFREQ)),0.,DLOG10(PHGEOSUM(IFREQ)))
            IF (PHGEOSEL(IFREQ).GT.0.D0)
     &        CALL hfillm
     &        (IDSEL+1,SNGL(FREQ(IFREQ)),0.,DLOG10(PHGEOSEL(IFREQ)))
            IF (PHBEAM(IFREQ).GT.0.D0)
     &        CALL hfillm
     &        (IDSEL+10,SNGL(FREQ(IFREQ)),0.,DLOG10(PHBEAM(IFREQ)))
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

          DO IFREQ=1,NFREQ
            CALL hfillm
     &        (IDSEL,SNGL(FREQ(IFREQ)),0.,PHGEOSUM(IFREQ))
            CALL hfillm
     &        (IDSEL+1,SNGL(FREQ(IFREQ)),0.,PHGEOSEL(IFREQ))
            CALL hfillm
     &        (IDSEL+10,SNGL(FREQ(IFREQ)),0.,PHBEAM(IFREQ))
          ENDDO

        ENDIF  !IFREQ2P

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
      DEALLOCATE(phspec)
      DEALLOCATE(phspecf)
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

      RETURN
      END
