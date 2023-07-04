*CMZ :  4.00/14 22/12/2021  18.07.25  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.15.16  by  Michael Scheer
*CMZ :  3.05/03 22/05/2018  07.13.27  by  Michael Scheer
*CMZ :  3.02/03 06/11/2014  14.54.13  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  11.18.26  by  Michael Scheer
*CMZ :  2.67/00 13/02/2012  10.58.17  by  Michael Scheer
*CMZ :  2.65/03 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.65/02 28/09/2009  13.02.09  by  Michael Scheer
*CMZ :  2.63/05 14/09/2009  15.19.42  by  Michael Scheer
*CMZ :  2.61/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  10.48.09  by  Michael Scheer
*CMZ :  2.51/01 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.50/00 29/04/2004  15.29.30  by  Michael Scheer
*CMZ :  2.41/10 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.37/07 11/12/2001  16.39.59  by  Michael Scheer
*CMZ :  2.37/01 14/11/2001  10.46.31  by  Michael Scheer
*CMZ :  2.37/00 13/11/2001  17.04.47  by  Michael Scheer
*CMZ :  2.36/00 08/11/2001  13.57.12  by  Michael Scheer
*CMZ :  2.34/09 21/09/2001  11.57.57  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  11.10.08  by  Michael Scheer
*CMZ :  2.33/00 03/05/2001  11.31.36  by  Michael Scheer
*CMZ :  2.31/01 25/04/2001  10.54.42  by  Michael Scheer
*CMZ :  2.31/00 24/04/2001  15.51.02  by  Michael Scheer
*CMZ :  2.30/03 20/04/2001  12.19.22  by  Michael Scheer
*CMZ :  2.30/02 12/04/2001  19.10.52  by  Michael Scheer
*CMZ :  2.30/01 12/04/2001  18.17.23  by  Michael Scheer
*CMZ :  2.20/12 11/04/2001  17.02.34  by  Michael Scheer
*CMZ :  2.20/11 11/04/2001  16.11.52  by  Michael Scheer
*CMZ :  2.20/10 10/04/2001  11.26.41  by  Michael Scheer
*CMZ :  2.20/09 03/04/2001  14.23.18  by  Michael Scheer
*CMZ :  2.15/01 30/03/2001  20.02.21  by  Michael Scheer
*CMZ :  2.20/08 18/03/2001  20.48.35  by  Michael Scheer
*CMZ :  2.20/07 18/03/2001  17.08.58  by  Michael Scheer
*CMZ :  2.20/06 15/03/2001  17.20.08  by  Michael Scheer
*CMZ :  2.20/05 15/03/2001  16.57.30  by  Michael Scheer
*CMZ :  2.20/04 09/03/2001  16.47.40  by  Michael Scheer
*CMZ :  2.20/03 23/02/2001  15.04.13  by  Michael Scheer
*CMZ :  2.20/02 21/02/2001  11.30.46  by  Michael Scheer
*CMZ :  2.20/01 20/02/2001  14.18.37  by  Michael Scheer
*-- Author : Michael Scheer

      SUBROUTINE SOUINTWIGS(ISOUR,IOBSV,INSIDE)
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

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEEP,workf90u.
      include 'workf90u.cmn'
*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEEP,amplif90u.
      include 'amplif90u.cmn'
*KEND.

C--- EVALUATE INTEGRALES FOR A SINGLE SOURCE

         IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,colli.
      include 'colli.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEND.

      COMPLEX*16 APOL
      COMPLEX*16 APOLH,APOLR,APOLL,APOL45

      INTEGER NFACP,IFAC,IFACO
      PARAMETER (NFACP=5)

      DOUBLE PRECISION AFREQR(3,NDFREQP),AFREQI(3,NDFREQP)
     &                    ,A1R(3),A1I(3),A2R(3),A2I(3)
     &  ,DGAMMA

        DOUBLE PRECISION DUMR,DUMI,FAC(NFACP)
      DOUBLE PRECISION AI1R,AE1R,AI2R,AE2R
        DOUBLE PRECISION EXPOMR,DEXPOMR,EXPOMI,DEXPOMI
      DOUBLE PRECISION AI1I,AE1I,AI2I,AE2I
     &                  ,BETP(3),BETP1(3),BETPP(3)

      DOUBLE PRECISION T0,T1,T2,TENDSOU,X0,X1,X2,X10,Y1,Y2,Z1,Z2,XENDSOU,R0
     &                  ,T,DT,DT2,DT0,DTIM00,DTIM01,VXP,VYP,VZP,TENDSOU1
     &                  ,R02,H2,H2R2,PHI,FREQR,CORRR0,R00
     &                  ,DMODUR,DMODUI,DMODU0R,DMODU0I
     &                  ,DDMODUR,DDMODUI
     &                  ,AXR,AXI,AYR,AYI,AZR,AZI
     &                  ,AX0R,AX0I,AY0R,AY0I,AZ0R,AZ0I
     &                  ,X2B,Y2B,Z2B
      DOUBLE PRECISION VX1,VY1,VZ1,BX1,BY1,BZ1
      DOUBLE PRECISION VX2,VY2,VZ2,BX2,BY2,BZ2,AX2D,AY2D,AZ2D
      DOUBLE PRECISION ECDUM,BS,BSQ,ECMAXS
      DOUBLE PRECISION TS,DPHASE,DPHSOUR(2,2)
        DOUBLE PRECISION GAMMA22,GAMMA221,GAMMA21,C1,OM,DOM,GAMMA,GAMGAM

        DOUBLE PRECISION BX,BY,BZ,RX,RY,RZ,PX,PY,PZ,RNBX,RNBY,RNBZ
      DOUBLE PRECISION R1,RNX,RNY,RNZ,DOM1,DOM2,BET1N,DUM11,R,BPX,BPY,BPZ
      DOUBLE PRECISION WGANG2,BET1N2,WGANG29
        DOUBLE PRECISION RARGOM,RARGDOM,RARG(5),PHASE,C
        double precision br2,rnr2,br4,rnr4,b3

      DOUBLE PRECISION, DIMENSION(3) :: RARGEXPOR,RARGEXPOI

        DOUBLE PRECISION EXPOMOR,EXPOMOI,DEXPOMOR,DEXPOMOI
      DOUBLE PRECISION EXPDOMR,EXPDOMI,DEXPDOMR,DEXPDOMI

      DOUBLE PRECISION DROIX,DTPHASE,DXEXI,CENXEXI

      INTEGER NOLDP
      PARAMETER (NOLDP=3)
        INTEGER IINSIDE,JINSIDE,IASYMP,IROIASY,INSIDE

        DOUBLE PRECISION STOK1,STOK2,STOK3,STOK4

      INTEGER ISOUR,IOBSV,IFREQ,IZAEHL,NZAEHL,IX10,I,ICAL,ICOMP
      INTEGER ICSPL

        INTEGER IROI,II,IZTOTS

      INTEGER NTUPP,IC,KWARN
        PARAMETER (NTUPP=29)
        REAL*8 FILLT(NTUPP)
        CHARACTER(4) CTUP(NTUPP)
        CHARACTER(4) ATUP(9)

        data ctup /'t','x','y','z','rx','ry','rz','rt','p','expr','expi','roi'
     &            ,'iob','ie','yob','zob','betn','dphi','dt','by2','isou'
     &            ,'spec','reax','imax','reay','imay','reaz','imaz','ifac'/

        data atup /'isou','iob','ie','icom','iae','a1r','a1i','a2r','a2i'/

      DATA ICAL/0/,KWARN/0/

        IF (IENELOSS.NE.0) STOP '*** SOUINTWIGS NOT READY FOR IENELOSS ***'




      IF (ICAL.EQ.0) THEN

       WRITE(LUNGFO,*)
       WRITE(LUNGFO,*)'       SUBROUTINE SOUINTWIGS:'
       WRITE(LUNGFO,*)

         IF (NFREQ.GT.NDFREQP) THEN
            WRITE(LUNGFO,*)
     &'*** ERROR IN SOUINTWIGS: NUMBER OF MAXIMUM PHOTON ENERGIES EXCEEDED'
            WRITE(LUNGFO,*)
     &'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
            WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED  ***'
            WRITE(6,*)
     &'*** ERROR IN SOUINTWIGS: NUMBER OF MAXIMUM PHOTON ENERGIES EXCEEDED'
            WRITE(6,*)
     &'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
            WRITE(6,*)'*** PROGRAM WAVE ABORTED  ***'
          STOP
         ENDIF    !(NFREQ.GT.NDFREQP)

          IROIASY=0
          DO IROI=1,NROIA
         IF (IROIASYEXP(IROI).NE.0) IROIASY=1
          ENDDO   !IROI

          IF (ISPECMODE.EQ.1) THEN
             DTIM00=DTMCO
         ELSE
             DTIM00=DTIM0
          ENDIF

          DTIM01=1.D0/DTIM00

            DO II=1,NSOURCE
            DO I=1,NROIA
                IWARNROI(I,II)=0
            ENDDO
            ENDDO

          GAMGAM=DMYGAMMA*DMYGAMMA
          GAMMA21=1.D0/GAMGAM
          GAMMA=DMYGAMMA
          GAMMA22=GAMGAM*2.D0
          GAMMA221=1.D0/GAMMA22
          WGANG2=(WGWINFC/GAMMA)**2+2.D0/(GAMGAM*(1.D0+DMYBETA))
          WGANG29=(WGWINFC/GAMMA*0.9)**2+2.D0/(GAMGAM*(1.D0+DMYBETA))
          DOM=(FREQ(2)-FREQ(1))/HBAREV1
          OM=FREQ(1)/HBAREV1
          C=CLIGHT1
          C1=1.D0/CLIGHT1

            IF (IWFILINT.LT.0) THEN
             CALL hbookm(NIDSOURCE,'RADIATION INTEGRALS',NTUPP
     &     , '//WAVE',nlpoi/jwfilint+2*jwfilint,CTUP)
           IF (IROIASY.NE.0) THEN
              CALL hbookm(NIDSOURCE+1,
     &       'ASYMPTOTIC EXPANSION OF RADIATION INTEGRALS',9
     &     ,  '//WAVE',nlpoi/jwfilint+2*jwfilint,ATUP)
           ENDIF   !(IROIASY.NE.0)
          ENDIF !(IWFILINT.LT.0)

          FAC(1)=1.D0/3.D0
          FAC(2)=4.D0/3.D0
          FAC(3)=2.D0/3.D0
          FAC(4)=4.D0/3.D0

          ICAL=1

      ENDIF !ICAL

      ICSPL=0
      FAC(5)=1.D0/3.D0

      IF (IOBSV.EQ.1) THEN
       WRITE(LUNGFO,*)'            SOURCE NUMBER',ISOUR,':'
       WRITE(LUNGFO,*)
      ENDIF

          IASYMP=0

          INSIDE=0
          IINSIDE=0
          JINSIDE=0
          SPECPOW(ISOUR+NSOURCE*(IOBSV-1))=0.D0

          EXPOMOR=1.D0
          EXPOMOI=0.D0
          DEXPOMOR=1.D0
          DEXPOMOI=0.D0

      DO IFREQ=1,NFREQ
          AFREQR(1,IFREQ)=0.D0
          AFREQR(2,IFREQ)=0.D0
          AFREQR(3,IFREQ)=0.D0
          AFREQI(1,IFREQ)=0.D0
          AFREQI(2,IFREQ)=0.D0
          AFREQI(3,IFREQ)=0.D0
      ENDDO   !IFREQ

      R0=OBSV(1,1)-SOURCEAO(1,1,ISOUR)

C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS     T=-R0*C1
        T=0.D0

      IF (ISPECMODE.EQ.1) THEN
          T0=DWT(1)
          T1=T0
          T2=DWT(MCO)
             XENDSOU=DWX(MCO)    !FINAL X
      ELSE
          T0=SOURCET(1,ISOUR)
          T1=T0
          T2=SOURCET(2,ISOUR)
             XENDSOU=SOURCEEO(1,1,ISOUR)    !FINAL X
      ENDIF

      TENDSOU=T2-T1

        X1=SOURCEAO(1,1,ISOUR)
        Y1=SOURCEAO(2,1,ISOUR)
        Z1=SOURCEAO(3,1,ISOUR)

        VX1=SOURCEAO(1,2,ISOUR)
        VY1=SOURCEAO(2,2,ISOUR)
        VZ1=SOURCEAO(3,2,ISOUR)

        BX1=SOURCEAO(1,4,ISOUR)
        BY1=SOURCEAO(2,4,ISOUR)
        BZ1=SOURCEAO(3,4,ISOUR)

      IF (IOBSV.EQ.ICBRILL) THEN
           ECSOUR(1,ISOUR)=0.0
           ECSOUR(2,ISOUR)=0.0
           ECSOUR(3,ISOUR)=0.0
           ECSOUR(4,ISOUR)=0.0
           ECMAX(   ISOUR)=-1.D30
      ENDIF

      ECMAXS=-1.D30
        IZTOTS=0

        X0=X1
        X2=X1
        X10=(XENDSOU-X0)/10.1D0

      NZAEHL=NLPOIO/NFACP*NFACP
      IF (NZAEHL.LT.NFACP) NZAEHL=NFACP
      NZAEHL=NZAEHL/2*2+1

      IF (ISPECMODE.EQ.1) THEN
          DT0=TENDSOU/(NZAEHL-1)
      ELSE
          DT0=TENDSOU/NZAEHL
      ENDIF

        DT=DT0

        IF (NROI.LT.0) THEN
            DROIX=(XENDSOU-X1)/(NROIA-1)
            DO IROI=1,NROIA
                ROIX(IROI)=X1+(IROI-1)*DROIX
                ROIP(IROI)=1.D0
            ENDDO
        ENDIF   !(NROI.LT.0)

        ROIX(1)=ROIX(1)-1.D-6
        ROIX(NROIA)=ROIX(NROIA)+1.D-6

        IF (X1.LT.ROIX(1).OR.XENDSOU.GT.ROIX(NROIA)) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN SOUINTWIGS: X OUTSIDE ROIS ***'
            WRITE(LUNGFO,*)'CHECK NAMELIST $ROIN'
            WRITE(LUNGFO,*)' *** PROGRAM WAVE ABORTED ***'
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN SOUINTWIGS: X OUTSIDE ROIS ***'
            WRITE(6,*)'CHECK NAMELIST $ROIN'
            WRITE(6,*)' *** PROGRAM WAVE ABORTED ***'
            STOP
        ENDIF   !IROI

      T=0.D0
      TS=0.D0

        X2=X1
        Y2=Y1
        Z2=Z1

        VX2=VX1
        VY2=VY1
        VZ2=VZ1

        BX2=BX1
        BY2=BY1
        BZ2=BZ1

C--- LOOP OVER STEPS

      PHASE=0.D0

        DO IROI=1,NROIA
            IPOIROI(IROI)=0
        ENDDO

        IROI=1
        DO I=1,NROIA
            IF (X1.GE.ROIX(I)) THEN
                IROI=I
            ENDIF !(X1.GE.ROIX(I))
        ENDDO   !IROI

      DT=DT0/ROIP(IROI)
      NZAEHL=NINT((TENDSOU-T)/DT)
      IF (NZAEHL.LT.NFACP) NZAEHL=NFACP
      NZAEHL=(NZAEHL+1)/NFACP*NFACP
      NZAEHL=NZAEHL/2*2+1
      IF (ISPECMODE.EQ.1) THEN
          DT=(TENDSOU-T)/(NZAEHL-1)
      ELSE
          DT=(TENDSOU-T)/NZAEHL
      ENDIF

      TENDSOU1=TENDSOU-DT
      DT2=DT/2.D0

C- CHECK STEPS SIZE

          IF (IWARNROI(IROI,ISOUR).EQ.0) THEN
          IF (DT.GT.DTIM00) THEN
                    WRITE(LUNGFO,*)
                    WRITE(LUNGFO,*)
     &'*** WARNING IN SOUINTWIGS, SOURCE, ROI:',ISOUR,IROI
                    WRITE(LUNGFO,*)
                    WRITE(LUNGFO,*)
     &             'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
                    WRITE(LUNGFO,*)'SIZE FOR TRAJECTORY!'
                    WRITE(LUNGFO,*)
                    WRITE(LUNGFO,*)
     &         'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
                    WRITE(6,*)
                    WRITE(6,*)
     &'*** WARNING IN SOUINTWIGS, SOURCE, ROI:',ISOUR,IROI
                    WRITE(6,*)
                    WRITE(6,*)'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
                    WRITE(6,*)'SIZE FOR TRAJECTORY!'
                    WRITE(6,*)
                    WRITE(6,*)
     &         'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
                    WRITE(6,*)
                   IWARNROI(IROI,ISOUR)=1
              ENDIF !DT
              ENDIF !IWARNROI

        IROI=IROI+1

        IX10=1
      IFAC=0

        IZAEHL=0 !LOOP COUNTER
      IF (IROIASY.EQ.0) THEN
          IZAEHL=IZAEHL+NOLDP
          NZAEHL=NZAEHL+NOLDP
      ENDIF

1000    IZAEHL=IZAEHL+1

        IF (IROI.LE.NROIA) THEN

            IF (X2.GE.ROIX(IROI)) THEN

            IF (IINSIDE.NE.2) THEN
            IF (IFACO.EQ.4) THEN

         FAC(5)=FAC(5)*(DT+DT0/ROIP(IROI))/DT

         IFAC=4

            ELSEIF (IFACO.EQ.5) THEN

C- CHECK STEPS SIZE
         IF (IROIASYEXP(IROI-1).NE.0) THEN
                    WRITE(LUNGFO,*)
                    WRITE(LUNGFO,*)
     &'*** WARNING IN SOUINTWIGS, SOURCE, ROI:',ISOUR,IROI
                    WRITE(LUNGFO,*)
                    WRITE(LUNGFO,*)
     &             'STEP SIZE CHANGED DURING ASYMPTOTIC EXPANSION PRECEDURE'
                    WRITE(LUNGFO,*)'THIS IS NOT RECOMMENDED, BE CAREFUL'
                    WRITE(LUNGFO,*)

                    WRITE(6,*)
                    WRITE(6,*)
     &'*** WARNING IN SOUINTWIGS, SOURCE, ROI:',ISOUR,IROI
                    WRITE(6,*)
                    WRITE(6,*)
     &             'STEP SIZE CHANGED DURING ASYMPTOTIC EXPANSION PRECEDURE'
                    WRITE(6,*)'THIS IS NOT RECOMMENDED, BE CAREFUL'
                    WRITE(6,*)
         ENDIF !IROIASYEXP


         DT=DT0/ROIP(IROI)
         NZAEHL=NINT((TENDSOU-T)/DT)
         IF (NZAEHL.LT.NFACP) NZAEHL=NFACP
         NZAEHL=(NZAEHL+1)/NFACP*NFACP
         NZAEHL=NZAEHL/2*2
         IF (ISPECMODE.EQ.1) THEN
             DT=(TENDSOU-T)/(NZAEHL-1)
         ELSE
             DT=(TENDSOU-T)/NZAEHL
         ENDIF

         TENDSOU1=TENDSOU-DT
         DT2=DT/2.D0

                 IF (IWARNROI(IROI,ISOUR).EQ.0) THEN

                 IF (DT.GT.DTIM00) THEN

                    WRITE(LUNGFO,*)
                    WRITE(LUNGFO,*)
     &'*** WARNING IN SOUINTWIGS, SOURCE, ROI:',ISOUR,IROI
                    WRITE(LUNGFO,*)
                    WRITE(LUNGFO,*)
     &             'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
                    WRITE(LUNGFO,*)'SIZE FOR TRAJECTORY!'
                    WRITE(LUNGFO,*)
                    WRITE(LUNGFO,*)
     &         'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
                    WRITE(6,*)
                    WRITE(6,*)
     &'*** WARNING IN SOUINTWIGS, SOURCE, ROI:',ISOUR,IROI
                    WRITE(6,*)
                    WRITE(6,*)'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
                    WRITE(6,*)'SIZE FOR TRAJECTORY!'
                    WRITE(6,*)
                    WRITE(6,*)
     &         'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
                    WRITE(6,*)

                   IWARNROI(IROI,ISOUR)=1

                ENDIF !DT
            ENDIF !IWARNROI

                IROI=IROI+1

         FAC(5)=FAC(1)
         IFAC=1

            ENDIF !IFAC
            ENDIF !IINSIDE

              ENDIF   !X2

        ENDIF   !IROI

        IPOIROI(IROI)=IPOIROI(IROI)+1

        T=T+DT

        X1=X2
        Y1=Y2
        Z1=Z2

        VX1=VX2
        VY1=VY2
        VZ1=VZ2

        BX1=BX2
        BY1=BY2
        BZ1=BZ2

      IF (ISPECMODE.NE.1) THEN

C GET MAGNETIC FIELD {

          X2B=X1+VX1*DT2
          Y2B=Y1+VY1*DT2
          Z2B=Z1+VZ1*DT2
          CALL MYBFELD(X2B,Y2B,Z2B,BX2,BY2,BZ2,AX2D,AY2D,AZ2D)

C GET MAGNETIC FIELD }

        BSQ=BX2*BX2+BY2*BY2+BZ2*BZ2
        BS=SQRT(BSQ)

C MOVE ONE STEP {

      CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DT,
     &             X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

          BX=VX2*C1
          BY=VY2*C1
          BZ=VZ2*C1

          BPX=VXP*C1
          BPY=VYP*C1
          BPZ=VZP*C1

C MOVE ONE STEP }

      ELSE  !ISPECMODE

          CALL WAVE_TRACK_INTER(TS,X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,BS,ICSPL,
     &      GAMMA)

      BSQ=BS*BS

      BX=VX2*C1
      BY=VY2*C1
      BZ=VZ2*C1
      BPX=VXP*C1
      BPY=VYP*C1
      BPZ=VZP*C1

      ENDIF !ISPECMODE

      BETP1(1)=BETP(1)
      BETP1(2)=BETP(2)
      BETP1(3)=BETP(3)
      BETP(1)=BPX
      BETP(2)=BPY
      BETP(3)=BPZ

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION {

C REAL PART OF INTEGRAND {

          RX=OBSV(1,IOBSV)-X2
          RY=OBSV(2,IOBSV)-Y2
          RZ=OBSV(3,IOBSV)-Z2

          R=SQRT(RX*RX+RY*RY+RZ*RZ)
          R1=1.D0/R

          RNX=RX*R1
          RNY=RY*R1
          RNZ=RZ*R1

C--- THE DISTANCE R IS INTRODUCED HERE EXPLICITLY (S. PROGRAM OF CHAOEN WANG

          BET1N=(1.D0-BX*RNX)-BY*RNY-BZ*RNZ
c 20090928{
      br2=by**2+bz**2
      rnr2=rny**2+rnz**2
      b3=dmybeta**3
      br4=br2**2
      rnr4=rnr2**2

      if(br2.lt.1.0d-4.and.rnr2.lt.1.0d-4) then
        bet1n=
     &    1.0d0/(1+dmybeta)/gamma**2
     &    +dmybeta*(rnr2/2.0d0
     &    +rnr4/8.0d0)
     &    +(br2/2.0d0
     &    -br2*rnr2/4.0d0
     &    -br2*rnr4/16.0d0)/dmybeta
     &    +b3*br4*(1.0d0/8.0d0
     &    -rnr2/16.0d0
     &    -rnr4/64.0d0)
     &    -by*rny
     &    -bz*rnz
      endif
c }20090928

          BET1N2=2.D0*BET1N
          DUM11=1.D0/BET1N
          DOM1=1.D0/(R*BET1N*BET1N)

          RNBX=RNX-BX
          RNBY=RNY-BY
          RNBZ=RNZ-BZ

          PX=(RNBY*BPZ-RNBZ*BPY)
          PY=(RNBZ*BPX-RNBX*BPZ)
          PZ=(RNBX*BPY-RNBY*BPX)

          IF (IVELOFIELD.EQ.0) THEN
              DOM2=C*DOM1*R1/GAMGAM
              RARG(1)=(RNY*PZ-RNZ*PY)*DOM1+(RNX-BX)*DOM2
              RARG(2)=(RNZ*PX-RNX*PZ)*DOM1+(RNY-BY)*DOM2
              RARG(3)=(RNX*PY-RNY*PX)*DOM1+(RNZ-BZ)*DOM2
            ELSE IF (IVELOFIELD.EQ.1) THEN
              RARG(1)=(RNY*PZ-RNZ*PY)*DOM1
              RARG(2)=(RNZ*PX-RNX*PZ)*DOM1
              RARG(3)=(RNX*PY-RNY*PX)*DOM1
            ELSE IF (IVELOFIELD.LT.0) THEN
              DOM2=C*DOM1*R1/GAMGAM
              RARG(1)=(RNX-BX)*DOM2
              RARG(2)=(RNY-BY)*DOM2
              RARG(3)=(RNZ-BZ)*DOM2
            ELSE  !IVELOFIELD
              WRITE(6,*)
     &          '*** ERROR IN SOUINTWIGS: BAD VALUE OF IVELOFIELD  ***'
              WRITE(6,*) '*** PROGRAM WAVE ABORTED  ***'
              STOP
            ENDIF !IVELOFIELD

          IF (IZAEHL.GE.NOLDP.AND.IINSIDE.EQ.0.AND.BET1N2.LE.WGANG2) THEN
         IF (IROIASYEXP(IROI-1).NE.0) IASYMP=1
         DPHSOUR(1,1)=BET1N*DT*FREQ(1)/HBAREV1
         DPHSOUR(1,2)=BET1N*DT*FREQ(NFREQ)/HBAREV1
         IINSIDE=1
         INSIDE=1
         JINSIDE=JINSIDE+1
         IF (JINSIDE.GT.1) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** WARNING IN SOUINTWIGS  ***'
           WRITE(LUNGFO,*)'*** SOURCE:',ISOUR
           WRITE(LUNGFO,*)'STRANGE SOURCE, CONTAINS SEVERAL SOURCES'
           WRITE(LUNGFO,*)'SOURCE AND OBSERVATION POINT:'
           WRITE(LUNGFO,*)ISOUR,OBSV(1,IOBSV),OBSV(2,IOBSV),OBSV(3,IOBSV)
           WRITE(LUNGFO,*)'AYMPTOTIC EXPANSION MAY BE INVALID!!'
           WRITE(LUNGFO,*)
     &           'RESULTS OF SPECTRUM CALCULATIONS MAY BE UNRELIABLE'
           WRITE(LUNGFO,*)'*** CHECK COLLIMATOR, PINHOLE, WGWINFC ... ***'
           WRITE(6,*)
           WRITE(6,*)'*** WARNING IN SOUINTWIGS  ***'
           WRITE(6,*)'*** SOURCE:',ISOUR
           WRITE(6,*)'*** STRANGE SOURCE, CONTAINS SEVERAL SOURCES'
           WRITE(6,*)'SOURCE AND OBSERVATION POINT:'
           WRITE(6,*)ISOUR,OBSV(1,IOBSV),OBSV(2,IOBSV),OBSV(3,IOBSV)
           WRITE(6,*)'*** CHECK COLLIMATOR, PINHOLE, WGWINFC ... ***'
           WRITE(6,*)'AYMPTOTIC EXPANSION IS INVALID!!'
           WRITE(6,*)'WARNING OF SPECTRUM CALCULATIONS ARE UNRELIABLE'
           JINSIDE=JINSIDE-1   !SUPRESS LOTS OF WARNINGS
           IASYMP=0  !BOUNDARY ALREADY TAKEN INTO ACCOUNT
         ENDIF !JINSIDE
              ELSEIF ((IINSIDE.EQ.1.AND.BET1N2.GT.WGANG2
     &          .OR.  IINSIDE.EQ.2).AND.IFAC.EQ.4) THEN
                IF (IROIASYEXP(IROI-1).NE.0) THEN
             IASYMP=2
         ELSEIF (IINSIDE.EQ.1) THEN
             IINSIDE=0
         ENDIF
         DPHSOUR(2,1)=BET1N*DT*FREQ(1)/HBAREV1
         DPHSOUR(2,2)=BET1N*DT*FREQ(NFREQ)/HBAREV1
          ELSE    !IINSIDE
         IASYMP=0
          ENDIF   !IINSIDE

          IF (IINSIDE.NE.0) THEN

C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS      RARG(4)=T+R*C1

C         RARG(4)=T-(RNX*X2+RNY*Y2+RNZ*Z2)*C1
C         RARG(4)=T+(R-R0)*C1

          DPHASE=BET1N*DT
          PHASE=PHASE+DPHASE
          RARG(4)=PHASE

          RARG(5)=
     &      (RARG(1)*RARG(1)+RARG(2)*RARG(2)+RARG(3)*RARG(3))*DUM11
C REAL PART OF INTEGRAND }

C COMPLEX PART OF INTEGRAND {

C    ASSUMES FREQ(I+1)=2*FREQ(I)   FOR IFREQ2P=2
C    OR FREQ(I+1)=FREQ(I)+DELTA    FOR IFREQ2P>2

C--- LOOP OVER ALL FREQUENCES

        IFREQ=1
        OM=FREQ(IFREQ)/HBAREV1

c       RARGOM=RARG(4)*OM
c       EXPOMR=COS(RARGOM)
c       EXPOMI=SIN(RARGOM)

        RARGOM=OM*DPHASE
        EXPDOMR=COS(RARGOM)
        EXPDOMI=SIN(RARGOM)
        DUMR=EXPOMOR*EXPDOMR-EXPOMOI*EXPDOMI
        DUMI=EXPOMOR*EXPDOMI+EXPOMOI*EXPDOMR
        EXPOMOR=DUMR
        EXPOMOI=DUMI
        EXPOMR=DUMR
        EXPOMI=DUMI

          IF(IFREQ2P.GT.2) THEN
c         RARGDOM=RARG(4)*DOM
c         DEXPOMR=COS(RARGDOM)
c         DEXPOMI=SIN(RARGDOM)
          RARGDOM=DOM*DPHASE
          DEXPDOMR=COS(RARGDOM)
          DEXPDOMI=SIN(RARGDOM)
          DUMR=DEXPOMOR*DEXPDOMR-DEXPOMOI*DEXPDOMI
          DUMI=DEXPOMOR*DEXPDOMI+DEXPOMOI*DEXPDOMR
          DEXPOMOR=DUMR
          DEXPOMOI=DUMI
          DEXPOMR=DUMR
          DEXPOMI=DUMI
        ENDIF  !IFREQ2P

          DO ICOMP=1,3
         RARGEXPOR(ICOMP)=RARG(ICOMP)*EXPOMR
          ENDDO
          DO ICOMP=1,3
         RARGEXPOI(ICOMP)=RARG(ICOMP)*EXPOMI
          ENDDO

         IF (IASYMP.EQ.1) THEN

          IF (BET1N2.LT.WGANG29) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** WARNING IN SOUINTWIGS  ***'
           WRITE(LUNGFO,*)'*** SOURCE:',ISOUR
           WRITE(LUNGFO,*)'STRANGE SOURCE, SMALL ENTRANCE ANGLE'
           WRITE(LUNGFO,*)'SOURCE AND OBSERVATION POINT:'
           WRITE(LUNGFO,*)ISOUR,IOBSV,SNGL(OBSV(1,IOBSV))
     &                          ,SNGL(OBSV(2,IOBSV)),SNGL(OBSV(3,IOBSV))
           WRITE(LUNGFO,*)'AYMPTOTIC EXPANSION MAY BE INVALID!!'
           WRITE(LUNGFO,*)
     &           'RESULTS OF SPECTRUM CALCULATIONS MAY BE UNRELIABLE'
           WRITE(LUNGFO,*)'*** CHECK COLLIMATOR, PINHOLE, WGWINFC ... ***'
           IF (KWARN.EQ.0) THEN
           WRITE(6,*)
           WRITE(6,*)'*** WARNING IN SOUINTWIGS  ***'
           WRITE(6,*)'STRANGE SOURCE, SMALL DEFLECTION ANGLE'
           WRITE(6,*)'AYMPTOTIC EXPANSION MIGHT BE INVALID!!'
           WRITE(6,*)'SEE OUTPUT FILE WAVE.OUT FOR DETAILS'
           KWARN=1
           ENDIF  !KWARN
          ENDIF   !WGANG29

             BETPP(1)=(BETP(1)-BETP1(1))/DT
             BETPP(2)=(BETP(2)-BETP1(2))/DT
             BETPP(3)=(BETP(3)-BETP1(3))/DT
                    CALL SOUASYEXP(IVELOFIELD,IROIASYEXP(IROI-1)
     &             ,X2,Y2,Z2,VX2,VY2,VZ2
     &             ,OBSV(1,IOBSV),OBSV(2,IOBSV),OBSV(3,IOBSV)
     &             ,BET1N,EXPOMR,EXPOMI,OM,BETP,BETPP,C,DMYGAMMA
     &             ,A1R,A1I,A2R,A2I)

                 DO ICOMP=1,3

             AI1R=A1R(ICOMP)
             AI1I=A1I(ICOMP)
             AI2R=A2R(ICOMP)
             AI2I=A2I(ICOMP)

             IF (IROIASYEXP(IROI-1).GT.0) THEN
              FILLT(1)=ISOUR
              FILLT(2)=IFREQ
              FILLT(3)=IOBSV
              FILLT(4)=ICOMP
                     FILLT(5)=1.0d0
              FILLT(6)=AI1R
              FILLT(7)=AI1I
              FILLT(8)=AI2R
              FILLT(9)=AI2I
              CALL hfm(NIDSOURCE+1,FILLT)
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
     &'     ind. of ener., obs. point , and comp.:'
     &                             ,IFREQ,IOBSV,ICOMP
              WRITE(LUNGFO,*)'     first order term at entrance:'
            WRITE(LUNGFO,*)'       ',AI1R,AI1I
              WRITE(LUNGFO,*)'     second order term at entrance:'
             WRITE(LUNGFO,*)'      ',AI2R,AI2I

             ENDIF   !(IROIASYEXP(IROI-1).GT.0)

             IF (X2.GE.XIANF) THEN
                  AFREQR(ICOMP,IFREQ)=AFREQR(ICOMP,IFREQ)+AI1R+AI2R
                  AFREQI(ICOMP,IFREQ)=AFREQI(ICOMP,IFREQ)+AI1I+AI2I
             ENDIF

           ENDDO  !ICOMP

         ELSEIF (IASYMP.EQ.2) THEN

          IF (BET1N2.LT.WGANG29) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** WARNING IN SOUINTWIGS  ***'
           WRITE(LUNGFO,*)'*** SOURCE:',ISOUR
           WRITE(LUNGFO,*)'STRANGE SOURCE, SMALL EXIT ANGLE'
           WRITE(LUNGFO,*)'SOURCE AND OBSERVATION POINT:'
           WRITE(LUNGFO,*)ISOUR,IOBSV,SNGL(OBSV(1,IOBSV))
     &                          ,SNGL(OBSV(2,IOBSV)),SNGL(OBSV(3,IOBSV))
           WRITE(LUNGFO,*)'AYMPTOTIC EXPANSION MAY BE INVALID!!'
           WRITE(LUNGFO,*)
     &           'RESULTS OF SPECTRUM CALCULATIONS MAY BE UNRELIABLE'
           WRITE(LUNGFO,*)'*** CHECK COLLIMATOR, PINHOLE, WGWINFC ... ***'
           IF (KWARN.EQ.0) THEN
           WRITE(6,*)
           WRITE(6,*)'*** WARNING IN SOUINTWIGS  ***'
           WRITE(6,*)'STRANGE SOURCE, SMALL DEFLECTION ANGLE'
           WRITE(6,*)'AYMPTOTIC EXPANSION MIGHT BE INVALID!!'
           WRITE(6,*)'SEE OUTPUT FILE WAVE.OUT FOR DETAILS'
           KWARN=1
           ENDIF  !KWARN
          ENDIF   !WGANG29

             BETPP(1)=(BETP(1)-BETP1(1))/DT
             BETPP(2)=(BETP(2)-BETP1(2))/DT
             BETPP(3)=(BETP(3)-BETP1(3))/DT
                    CALL SOUASYEXP(IVELOFIELD,IROIASYEXP(IROI-1)
     &             ,X2,Y2,Z2,VX2,VY2,VZ2
     &             ,OBSV(1,IOBSV),OBSV(2,IOBSV),OBSV(3,IOBSV)
     &             ,BET1N,EXPOMR,EXPOMI,OM,BETP,BETPP,C,DMYGAMMA
     &             ,A1R,A1I,A2R,A2I)

                 DO ICOMP=1,3

             AE1R=A1R(ICOMP)
             AE1I=A1I(ICOMP)
             AE2R=A2R(ICOMP)
             AE2I=A2I(ICOMP)

             IF (X2.GE.XIANF) THEN
                  AFREQR(ICOMP,IFREQ)=AFREQR(ICOMP,IFREQ)-AE1R-AE2R
                  AFREQI(ICOMP,IFREQ)=AFREQI(ICOMP,IFREQ)-AE1I-AE2I
                  ENDIF

             IF (IROIASYEXP(IROI-1).GT.0) THEN
              FILLT(1)=ISOUR
              FILLT(2)=IFREQ
              FILLT(3)=IOBSV
              FILLT(4)=ICOMP
                     FILLT(5)=2.0d0
              FILLT(6)=AE1R
              FILLT(7)=AE1I
              FILLT(8)=AE2R
              FILLT(9)=AE2I
              CALL hfm(NIDSOURCE+1,FILLT)

              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
     &'     ind. of ener., obs. point , and comp.:'
     &                             ,IFREQ,IOBSV,ICOMP
              WRITE(LUNGFO,*)'     first order term at exit:'
            WRITE(LUNGFO,*)'      ',AE1R,AE1I
              WRITE(LUNGFO,*)'     second order term at exit:'
            WRITE(LUNGFO,*)'      ',AE2R,AE2I
             ENDIF   !(IROIASYEXP(IROI-1).GT.0)

           ENDDO  !ICOMP

         ENDIF !IASYMP

        IF (X2.GE.XIANF) THEN

          IF (IOBSV.EQ.ICBRILL) THEN
               ECSOUR(1,ISOUR)=ECSOUR(1,ISOUR)+BS
               ECSOUR(4,ISOUR)=ECSOUR(4,ISOUR)+SIGN(BS,BY2)
               ECSOUR(3,ISOUR)=ECSOUR(3,ISOUR)+BSQ
               IF (ECMAX(ISOUR).LT.BS) ECMAX(ISOUR)=BS
               IZTOTS=IZTOTS+1
          ENDIF
             IF (ECMAXS.LT.BS) ECMAXS=BS

            IFAC=IFAC+1
            ILIOB=ISOUR+NSOURCE*(IOBSV-1)
              SPECPOW(ILIOB)=SPECPOW(ILIOB)+FAC(IFAC)*RARG(5)*DT
               DO ICOMP=1,3
                 AFREQR(ICOMP,IFREQ)=AFREQR(ICOMP,IFREQ)
     &         +FAC(IFAC)*RARGEXPOR(ICOMP)*DT
               ENDDO   !ICOMP
               DO ICOMP=1,3
                 AFREQI(ICOMP,IFREQ)=AFREQI(ICOMP,IFREQ)
     &         +FAC(IFAC)*RARGEXPOI(ICOMP)*DT
               ENDDO   !ICOMP
        ENDIF  !XIANF

        IF (IWFILINT.NE.0) THEN
        IF (IWFILINT.LT.0) THEN
            FILLT(1)=T
            FILLT(2)=X2
            FILLT(3)=Y2
            FILLT(4)=Z2
            FILLT(5)=RARG(1)
            FILLT(6)=RARG(2)
            FILLT(7)=RARG(3)
            FILLT(8)=RARG(4)
            FILLT(9)=RARG(5)
            FILLT(10)=EXPOMR
            FILLT(11)=EXPOMI
            FILLT(12)=IROI-1
            FILLT(13)=IOBSV
            FILLT(14)=IFREQ
            FILLT(15)=OBSV(2,IOBSV)
            FILLT(16)=OBSV(3,IOBSV)
            FILLT(17)=BET1N
            FILLT(18)=DPHASE*OM
            FILLT(19)=DT
            FILLT(20)=BY2
            FILLT(21)=ISOUR
            FILLT(22)=
     &         (
     &          AFREQR(1,IFREQ)*AFREQR(1,IFREQ)
     &         +AFREQI(1,IFREQ)*AFREQI(1,IFREQ)
     &         +AFREQR(2,IFREQ)*AFREQR(2,IFREQ)
     &         +AFREQI(2,IFREQ)*AFREQI(2,IFREQ)
     &         +AFREQR(3,IFREQ)*AFREQR(3,IFREQ)
     &         +AFREQI(3,IFREQ)*AFREQI(3,IFREQ)
     &         )*SPECNOR
            FILLT(23)=AFREQR(1,IFREQ)*SPECNOR
            FILLT(24)=AFREQI(1,IFREQ)*SPECNOR
            FILLT(25)=AFREQR(2,IFREQ)*SPECNOR
            FILLT(26)=AFREQI(2,IFREQ)*SPECNOR
            FILLT(27)=AFREQR(3,IFREQ)*SPECNOR
            FILLT(28)=AFREQI(3,IFREQ)*SPECNOR
            FILLT(29)=IFAC

            CALL hfm(NIDSOURCE,FILLT)

        ELSE IF (ISOUR.EQ.IWFILINT.AND.IOBSV.EQ.1) THEN

         WRITE(LUNINT,*) IZAEHL,IFREQ,X2
         WRITE(LUNINT,*) (RARG(1),IC=1,3)
         WRITE(LUNINT,*) RARG(4)*OM,RARG(5)
         WRITE(LUNINT,*)EXPOMR,EXPOMI
         WRITE(LUNINT,*)RARG(1)*EXPOMR,RARG(1)*EXPOMI
         WRITE(LUNINT,*)RARG(2)*EXPOMR,RARG(2)*EXPOMI
         WRITE(LUNINT,*)RARG(3)*EXPOMR,RARG(3)*EXPOMI

        ENDIF !IWFILINT.LT.0
        ENDIF !IWFILINT.NE.0

          DO IFREQ=2,NFREQ

          OM=OM+DOM

          IF    (IFREQ2P.GT.2) THEN
         DUMR=EXPOMR*DEXPOMR-EXPOMI*DEXPOMI
         DUMI=EXPOMR*DEXPOMI+EXPOMI*DEXPOMR
         EXPOMR=DUMR
         EXPOMI=DUMI
          ELSEIF(IFREQ2P.EQ.2) THEN
         DUMR=EXPOMR*EXPOMR-EXPOMI*EXPOMI
         DUMI=EXPOMR*EXPOMI+EXPOMI*EXPOMR
         EXPOMR=DUMR
         EXPOMI=DUMI
          ELSE
         OM=FREQ(IFREQ)/HBAREV1
         EXPOMR=COS(RARG(4)*OM)
         EXPOMI=SIN(RARG(4)*OM)
          ENDIF

           DO ICOMP=1,3
          RARGEXPOR(ICOMP)=RARG(ICOMP)*EXPOMR
           ENDDO
           DO ICOMP=1,3
          RARGEXPOI(ICOMP)=RARG(ICOMP)*EXPOMI
           ENDDO

         IF (IASYMP.EQ.1) THEN

             BETPP(1)=(BETP(1)-BETP1(1))/DT
             BETPP(2)=(BETP(2)-BETP1(2))/DT
             BETPP(3)=(BETP(3)-BETP1(3))/DT
                    CALL SOUASYEXP(IVELOFIELD,IROIASYEXP(IROI-1)
     &             ,X2,Y2,Z2,VX2,VY2,VZ2
     &             ,OBSV(1,IOBSV),OBSV(2,IOBSV),OBSV(3,IOBSV)
     &             ,BET1N,EXPOMR,EXPOMI,OM,BETP,BETPP,C,DMYGAMMA
     &             ,A1R,A1I,A2R,A2I)

              DO ICOMP=1,3

             AI1R=A1R(ICOMP)
             AI1I=A1I(ICOMP)
             AI2R=A2R(ICOMP)
             AI2I=A2I(ICOMP)

             IF (IROIASYEXP(IROI-1).GT.0) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
     &'     ind. of ener., obs. point , and comp.:'
     &                             ,IFREQ,IOBSV,ICOMP
              WRITE(LUNGFO,*)'     first order term at entrance:'
            WRITE(LUNGFO,*)'      ',AI1R,AI1I
              WRITE(LUNGFO,*)'     second order term at entrance:'
            WRITE(LUNGFO,*)'      ',AI2R,AI2I
             ENDIF   !(IROIASYEXP(IROI-1).GT.0)

             IF (X2.GE.XIANF) THEN
                  AFREQR(ICOMP,IFREQ)=AFREQR(ICOMP,IFREQ)+AI1R+AI2R
                  AFREQI(ICOMP,IFREQ)=AFREQI(ICOMP,IFREQ)+AI1I+AI2I
             ENDIF

          ENDDO   !ICOMP

         ELSEIF (IASYMP.EQ.2) THEN

             BETPP(1)=(BETP(1)-BETP1(1))/DT
             BETPP(2)=(BETP(2)-BETP1(2))/DT
             BETPP(3)=(BETP(3)-BETP1(3))/DT
                    CALL SOUASYEXP(IVELOFIELD,IROIASYEXP(IROI-1)
     &             ,X2,Y2,Z2,VX2,VY2,VZ2
     &             ,OBSV(1,IOBSV),OBSV(2,IOBSV),OBSV(3,IOBSV)
     &             ,BET1N,EXPOMR,EXPOMI,OM,BETP,BETPP,C,DMYGAMMA
     &             ,A1R,A1I,A2R,A2I)

              DO ICOMP=1,3

             AE1R=A1R(ICOMP)
             AE1I=A1I(ICOMP)
             AE2R=A2R(ICOMP)
             AE2I=A2I(ICOMP)

             IF (IROIASYEXP(IROI-1).GT.0) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
     &'     ind. of ener., obs. point , and comp.:'
     &                             ,IFREQ,IOBSV,ICOMP
              WRITE(LUNGFO,*)'     first order term at exit:'
            WRITE(LUNGFO,*)'      ',AE1R,AE1I
              WRITE(LUNGFO,*)'     second order term exit:'
             WRITE(LUNGFO,*)'      ',AE2R,AE2I
             ENDIF   !(IROIASYEXP(IROI-1).GT.0)

             IF (X2.GE.XIANF) THEN
                  AFREQR(ICOMP,IFREQ)=AFREQR(ICOMP,IFREQ)-AE1R-AE2R
                  AFREQI(ICOMP,IFREQ)=AFREQI(ICOMP,IFREQ)-AE1I-AE2I
             ENDIF

           ENDDO  !ICOMP

         ENDIF !(IASYMP.EQ.2)

        IF (X2.GE.XIANF) THEN
               DO ICOMP=1,3
                 AFREQR(ICOMP,IFREQ)=AFREQR(ICOMP,IFREQ)
     &         +FAC(IFAC)*RARGEXPOR(ICOMP)*DT
               ENDDO   !ICOMP
               DO ICOMP=1,3
                 AFREQI(ICOMP,IFREQ)=AFREQI(ICOMP,IFREQ)
     &         +FAC(IFAC)*RARGEXPOI(ICOMP)*DT
               ENDDO   !ICOMP
        ENDIF

        IF (IWFILINT.NE.0) THEN
        IF (IWFILINT.LT.0) THEN
            FILLT(1)=T
            FILLT(2)=X2
            FILLT(3)=Y2
            FILLT(4)=Z2
            FILLT(5)=RARG(1)
            FILLT(6)=RARG(2)
            FILLT(7)=RARG(3)
            FILLT(8)=RARG(4)
            FILLT(9)=RARG(5)
            FILLT(10)=EXPOMR
            FILLT(11)=EXPOMI
            FILLT(12)=IROI-1
            FILLT(13)=IOBSV
            FILLT(14)=IFREQ
            FILLT(15)=OBSV(2,IOBSV)
            FILLT(16)=OBSV(3,IOBSV)
            FILLT(17)=BET1N
            FILLT(18)=DPHASE*OM
            FILLT(19)=DT
            FILLT(20)=BY2
            FILLT(21)=ISOUR
            FILLT(22)=
     &         (
     &          AFREQR(1,IFREQ)*AFREQR(1,IFREQ)
     &         +AFREQI(1,IFREQ)*AFREQI(1,IFREQ)
     &         +AFREQR(2,IFREQ)*AFREQR(2,IFREQ)
     &         +AFREQI(2,IFREQ)*AFREQI(2,IFREQ)
     &         +AFREQR(3,IFREQ)*AFREQR(3,IFREQ)
     &         +AFREQI(3,IFREQ)*AFREQI(3,IFREQ)
     &         )*SPECNOR
            FILLT(23)=AFREQR(1,IFREQ)*SPECNOR
            FILLT(24)=AFREQI(1,IFREQ)*SPECNOR
            FILLT(25)=AFREQR(2,IFREQ)*SPECNOR
            FILLT(26)=AFREQI(2,IFREQ)*SPECNOR
            FILLT(27)=AFREQR(3,IFREQ)*SPECNOR
            FILLT(28)=AFREQI(3,IFREQ)*SPECNOR
            FILLT(29)=IFAC
            CALL hfm(NIDSOURCE,FILLT)

        ELSEIF (ISOUR.EQ.IWFILINT.AND.IOBSV.EQ.1) THEN

         WRITE(LUNINT,*) IZAEHL,IFREQ,X2
         WRITE(LUNINT,*) (RARG(1),IC=1,3)
         WRITE(LUNINT,*) RARG(4)*OM,RARG(5)
         WRITE(LUNINT,*)EXPOMR,EXPOMI
         WRITE(LUNINT,*)RARG(1)*EXPOMR,RARG(1)*EXPOMI
         WRITE(LUNINT,*)RARG(2)*EXPOMR,RARG(2)*EXPOMI
         WRITE(LUNINT,*)RARG(3)*EXPOMR,RARG(3)*EXPOMI

        ENDIF !IWFILINT.LT.0
        ENDIF !IWFILINT.NE.0

        ENDDO   !LOOP OVER ALL FREQUENCES

      IFACO=IFAC
      IF (IFAC.EQ.4) IFAC=2

        ENDIF   !IINSIDE

C COMPLEX PART OF INTEGRAND }

      IF (IASYMP.EQ.2) THEN
          IINSIDE=0
      ENDIF

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION }

        TS=TS+DT

C--- END OF LOOP OVER TIME STEPS

      IF (T.LT.TENDSOU1.AND.X2.LE.XIEND)  GOTO 1000

      IF (IFAC.NE.5.AND.IFAC.NE.0) THEN
          IF (IFACO.EQ.4) IFAC=4
        IINSIDE=2
        GOTO 1000
      ENDIF

      IF (IINSIDE.NE.0) THEN
         DPHSOUR(2,1)=BET1N*DT*FREQ(1)/HBAREV1
         DPHSOUR(2,2)=BET1N*DT*FREQ(NFREQ)/HBAREV1
      ENDIF

C- STORE NUMBER OF POINTS FOR INTEGRATION

        IF (IOBSV.EQ.ICBRILL) IPOISOU(ISOUR)=IZAEHL

      IF (IOBSV.EQ.ICBRILL.AND.IZTOTS.GT.0) THEN

        ECSOUR(1,ISOUR)=ECSOUR(1,ISOUR)/IZTOTS
        ECSOUR(4,ISOUR)=ECSOUR(4,ISOUR)/IZTOTS
        ECSOUR(3,ISOUR)=ECSOUR(3,ISOUR)/IZTOTS

        ECDUM=ECSOUR(3,ISOUR)-ECSOUR(1,ISOUR)**2

        IF (ECDUM.LT.0.0) ECDUM=0.

      IF (ECSOUR(1,ISOUR).NE.0.D0) THEN
           ECSOUR(3,ISOUR)=SQRT(ECDUM)/ECSOUR(1,ISOUR)
      ELSE
        ECSOUR(3,ISOUR)=0.D0
      ENDIF

        ECSOUR(2,ISOUR)=ECSOUR(1,ISOUR)*ecdipev1*DMYENERGY**2   !CRITICAL ENERGY

      ENDIF !ISPECMODE

      IF (IAMPLI.LT.0) THEN
         DXEXI=MIN(SOURCEEO(1,1,ISOUR),XIEND)
     &               -MAX(SOURCEAO(1,1,ISOUR),XIANF)
         CENXEXI=(MIN(SOURCEEO(1,1,ISOUR),XIEND)
     &                 +MAX(SOURCEAO(1,1,ISOUR),XIANF))/2.D0
         DTPHASE=(WTRA2IS(ISOUR)+GAMMA21*DXEXI/2.D0)/CLIGHT1
            FREQR=2.D0*PI1/DTPHASE*HBAREV1
      ENDIF !(IAMPLI.LT.0)

      DO IFREQ=1,NFREQ

        ILIOBFR=ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1))
          IFROB=IFREQ+NFREQ*(IOBSV-1)
          IOBFR=IOBSV+NOBSV*(IFREQ-1)

        IF (IAMPLI.LT.0) THEN

          AX0R=AFREQR(1,IFREQ)
          AX0I=AFREQI(1,IFREQ)
          AY0R=AFREQR(2,IFREQ)
          AY0I=AFREQI(2,IFREQ)
          AZ0R=AFREQR(3,IFREQ)
          AZ0I=AFREQI(3,IFREQ)

          AXR=AX0R
          AXI=AX0I
          AYR=AY0R
          AYI=AY0I
          AZR=AZ0R
          AZI=AZ0I

          AFREQR(1,IFREQ)=0.D0
          AFREQI(1,IFREQ)=0.D0
          AFREQR(2,IFREQ)=0.D0
          AFREQI(2,IFREQ)=0.D0
          AFREQR(3,IFREQ)=0.D0
          AFREQI(3,IFREQ)=0.D0

            R0=OBSV(1,NOBSV/2+1)-CENXEXI
            R02=R0*R0
          R00=R0
            H2=(OBSV(2,IOBSV))**2+(OBSV(3,IOBSV))**2
            H2R2=H2/R02

          DTPHASE=(WTRA2IS(ISOUR)+(H2R2+GAMMA21)*DXEXI/2.D0)/CLIGHT1
          PHI=2.D0*PI1*FREQ(IFREQ)*ECHARGE1/HPLANCK1*DTPHASE

          DMODUR=COS(PHI)
          DMODUI=SIN(PHI)
            DMODU0R=DMODUR
            DMODU0I=DMODUI
          DDMODUR=1.D0
          DDMODUI=0.D0

          DO I=1,-IAMPLI

          R0=OBSV(1,NOBSV/2+1)+DXEXI/2.D0*(-IAMPLI-2*(I-1))
          CORRR0=R00/R0
            R02=R0*R0
            H2=(OBSV(2,IOBSV))**2+(OBSV(3,IOBSV))**2
            H2R2=H2/R02

          DTPHASE=(WTRA2IS(ISOUR)+(H2R2+GAMMA21)*DXEXI/2.D0)/CLIGHT1
          PHI=2.D0*PI1*FREQ(IFREQ)*ECHARGE1/HPLANCK1*DTPHASE

          DMODUR=COS(PHI)
          DMODUI=SIN(PHI)
            DMODU0R=DMODUR
            DMODU0I=DMODUI
          DDMODUR=1.D0
          DDMODUI=0.D0

          AFREQR(1,IFREQ)=AFREQR(1,IFREQ)+AXR
          AFREQI(1,IFREQ)=AFREQI(1,IFREQ)+AXI
          AFREQR(2,IFREQ)=AFREQR(2,IFREQ)+AYR
          AFREQI(2,IFREQ)=AFREQI(2,IFREQ)+AYI
          AFREQR(3,IFREQ)=AFREQR(3,IFREQ)+AZR
          AFREQI(3,IFREQ)=AFREQI(3,IFREQ)+AZI

          IF (AMPRAN.NE.0.D0) THEN

             PHI=2.D0*PI1*XRANA(I)/FREQR*FREQ(IFREQ)

             DDMODUR=COS(PHI)
             DDMODUI=SIN(PHI)

          ENDIF   !(AMPRAN.NE.0.D0)

          DUMR=AX0R*DMODU0R-AX0I*DMODU0I
          DUMI=AX0R*DMODU0I+AX0I*DMODU0R
          AX0R=DUMR
          AX0I=DUMI
          DUMR=AY0R*DMODU0R-AY0I*DMODU0I
          DUMI=AY0R*DMODU0I+AY0I*DMODU0R
          AY0R=DUMR
          AY0I=DUMI
          DUMR=AZ0R*DMODU0R-AZ0I*DMODU0I
          DUMI=AZ0R*DMODU0I+AZ0I*DMODU0R
          AZ0R=DUMR
          AZ0I=DUMI
          AXR=AX0R*CORRR0
          AXI=AX0I*CORRR0
          AYR=AY0R*CORRR0
          AYI=AY0I*CORRR0
          AZR=AZ0R*CORRR0
          AZI=AZ0I*CORRR0

          DMODUR=DMODU0R*DDMODUR-DMODU0I*DDMODUI
          DMODUI=DMODU0R*DDMODUI+DMODU0I*DDMODUR

          DUMR=AXR*DMODUR-AXI*DMODUI
          DUMI=AXR*DMODUI+AXI*DMODUR
          AXR=DUMR
          AXI=DUMI
          DUMR=AYR*DMODUR-AYI*DMODUI
          DUMI=AYR*DMODUI+AYI*DMODUR
          AYR=DUMR
          AYI=DUMI
          DUMR=AZR*DMODUR-AZI*DMODUI
          DUMI=AZR*DMODUI+AZI*DMODUR
          AZR=DUMR
          AZI=DUMI

        ENDDO !IAMPLI

        ENDIF  !(IAMPLI.LT.0)

        AFREQ(1,IFROB)=DCMPLX(AFREQR(1,IFREQ),AFREQI(1,IFREQ))*REFLEC(1)
        AFREQ(2,IFROB)=DCMPLX(AFREQR(2,IFREQ),AFREQI(2,IFREQ))*REFLEC(2)
        AFREQ(3,IFROB)=DCMPLX(AFREQR(3,IFREQ),AFREQI(3,IFREQ))*REFLEC(3)


        IF(SPECCUT.GT.0.0D0) THEN
        IF (ISPECMODE.EQ.1) ECMAXS=ECMAX(ISOUR)
        IF(FREQ(IFREQ).GT.SPECCUT*ecdipev1*DMYENERGY**2*ECMAXS) THEN
          AFREQ(1,IFROB)=(0.D0,0.D0)
          AFREQ(2,IFROB)=(0.D0,0.D0)
          AFREQ(3,IFROB)=(0.D0,0.D0)
        ENDIF
        ENDIF

          IF (IPOLA.EQ.0) THEN

            SPEC(ILIOBFR)=
     &        DREAL(
     &          AFREQ(1,IFROB)*CONJG(AFREQ(1,IFROB))
     &         +AFREQ(2,IFROB)*CONJG(AFREQ(2,IFROB))
     &         +AFREQ(3,IFROB)*CONJG(AFREQ(3,IFROB))
     &         )*SPECNOR

            REAIMA(1,1,IOBFR)=REAIMA(1,1,IOBFR)+
     &        DREAL(AFREQ(1,IFROB))
            REAIMA(2,1,IOBFR)=REAIMA(2,1,IOBFR)+
     &        DREAL(AFREQ(2,IFROB))
            REAIMA(3,1,IOBFR)=REAIMA(3,1,IOBFR)+
     &        DREAL(AFREQ(3,IFROB))

            REAIMA(1,2,IOBFR)=REAIMA(1,2,IOBFR)+
     &        DIMAG(AFREQ(1,IFROB))
            REAIMA(2,2,IOBFR)=REAIMA(2,2,IOBFR)+
     &        DIMAG(AFREQ(2,IFROB))
            REAIMA(3,2,IOBFR)=REAIMA(3,2,IOBFR)+
     &        DIMAG(AFREQ(3,IFROB))

          ELSE    !IPOLA

          APOL=
     &          AFREQ(1,IFROB)*CONJG(VPOLA(1))
     &         +AFREQ(2,IFROB)*CONJG(VPOLA(2))
     &         +AFREQ(3,IFROB)*CONJG(VPOLA(3))

            SPEC(ILIOBFR)=
     &        DREAL(APOL*CONJG(APOL))*SPECNOR
            REAIMA(1,1,IOBFR)=REAIMA(1,1,IOBFR)+
     &        DREAL(AFREQ(1,IFROB))
            REAIMA(2,1,IOBFR)=REAIMA(2,1,IOBFR)+
     &        DREAL(AFREQ(2,IFROB))
            REAIMA(3,1,IOBFR)=REAIMA(3,1,IOBFR)+
     &        DREAL(AFREQ(3,IFROB))

            REAIMA(1,2,IOBFR)=REAIMA(1,2,IOBFR)+
     &        DIMAG(AFREQ(1,IFROB))
            REAIMA(2,2,IOBFR)=REAIMA(2,2,IOBFR)+
     &        DIMAG(AFREQ(2,IFROB))
            REAIMA(3,2,IOBFR)=REAIMA(3,2,IOBFR)+
     &        DIMAG(AFREQ(3,IFROB))

          ENDIF   !IPOLA

      IF (ISTOKES.NE.0) THEN

          APOLH=
     &          AFREQ(1,IFROB)*CONJG(VSTOKES(1,1))
     &         +AFREQ(2,IFROB)*CONJG(VSTOKES(1,2))
     &         +AFREQ(3,IFROB)*CONJG(VSTOKES(1,3))

          APOLR=
     &          AFREQ(1,IFROB)*CONJG(VSTOKES(2,1))
     &         +AFREQ(2,IFROB)*CONJG(VSTOKES(2,2))
     &         +AFREQ(3,IFROB)*CONJG(VSTOKES(2,3))

          APOLL=
     &          AFREQ(1,IFROB)*CONJG(VSTOKES(3,1))
     &         +AFREQ(2,IFROB)*CONJG(VSTOKES(3,2))
     &         +AFREQ(3,IFROB)*CONJG(VSTOKES(3,3))

          APOL45=
     &          AFREQ(1,IFROB)*CONJG(VSTOKES(4,1))
     &         +AFREQ(2,IFROB)*CONJG(VSTOKES(4,2))
     &         +AFREQ(3,IFROB)*CONJG(VSTOKES(4,3))

            STOK1=
     &        APOLR*CONJG(APOLR)+
     &        APOLL*CONJG(APOLL)

            STOK2=-STOK1+
     &        2.*APOLH*CONJG(APOLH)

            STOK3=
     &        2.*APOL45*CONJG(APOL45)-
     &        STOK1

            STOK4=
     &        APOLR*CONJG(APOLR)-
     &        APOLL*CONJG(APOLL)


            STOKES(1,IOBFR)=STOKES(1,IOBFR)+
     &                          STOK1*SPECNOR

            STOKES(2,IOBFR)=STOKES(2,IOBFR)+
     &                          STOK2*SPECNOR

            STOKES(3,IOBFR)=STOKES(3,IOBFR)+
     &                          STOK3*SPECNOR

            STOKES(4,IOBFR)=STOKES(4,IOBFR)+
     &                          STOK4*SPECNOR

      ENDIF !ISTOKES

      ENDDO !IFREQ

        ILIOB=ISOUR+NSOURCE*(IOBSV-1)
        SPECPOW(ILIOB)=SPECPOW(ILIOB)
     &                *ECHARGE1/16.D0/PI1/PI1/EPS01/C
     &                *DMYCUR     !NUMBER OF e-

      IF (IOBSV.EQ.ICBRILL) THEN

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &'       phase advance per step at beginning and end of source for'
        WRITE(LUNGFO,*)
     &'       lowest and highest photon energy at selected observation point:'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'       beginning:',SNGL(DPHSOUR(1,1)),SNGL(DPHSOUR(1,2))
        WRITE(LUNGFO,*)'       end:      ',SNGL(DPHSOUR(2,1)),SNGL(DPHSOUR(2,2))
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'       ROIs (boundary, precision, points, asymp. exp.:'
        WRITE(LUNGFO,*)

        DO IROI=1,NROIA-1
                 WRITE(LUNGFO,*)
     &           IROI,SNGL(ROIX(IROI)),SNGL(ROIP(IROI)),IPOIROI(IROI+1)
     &          ,IROIASYEXP(IROI)
        ENDDO
                 WRITE(LUNGFO,*)
     &           NROI,SNGL(ROIX(NROIA))

      ENDIF !IOBSV

        IF (IOBSV.EQ.NOBSV) THEN

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'       SOURCE, TOTAL NUMBER OF STEPS:',ISOUR,IZAEHL
        WRITE(LUNGFO,*)'       (controlled by NLPOI and namelist $ROIN)'
        WRITE(LUNGFO,*)

      ENDIF

      RETURN
      END
