*CMZ :  3.08/01 02/04/2019  15.33.15  by  Michael Scheer
*CMZ :  3.02/00 18/09/2014  14.21.18  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE SPECDIP4
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
*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

C--- FILL ARRAY SPEC WITH DIPOL SPECTRUM FUNCTION

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,colli.
      include 'colli.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,specdip.
      include 'specdip.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,datetime.
      include 'datetime.cmn'
*KEND.

      INTEGER ISOUR,IFREQ,IOBSV,ISPLN,ITANG,I,IC,JC,IFAIL,J,JDX10,JX10,IX10
     &  ,IWARN,IEPS,IWARNOB

      DOUBLE PRECISION T,XT,YT,ZT,VXT,VYT,VZT,VXP,VYP,VZP,BS
      DOUBLE PRECISION RX,RY,RZ,RR3,RR,PSI,DFDTDP,PAR,PER,POWR,OBANG,OBANGMN
     &  ,ANG3(3),XT3(3),YT3(3),ZT3(3),VXT3(3),VYT3(3),VZT3(3),A(3),
     &  dum3(3)
      DOUBLE PRECISION Y,OMEGAC,BX,BY,BZ,BRX,BRY,BRZ,RN,VN,DUM,BRN

      DOUBLE PRECISION DX2,DZY2,EPS(6),DRRED,ANS,OM,
     &  APERANG,APERV(3),APERCORR,
     &  RXX,RYY,RZZ,RRN,RRR,
     &  GAMMA,speck

      REAL STOK1,STOK2,STOK3,STOK4

      COMPLEX*16 APOL,AFREQ(3),EXPOM
      COMPLEX*8 APOLH,APOLR,APOLL,APOL45

      DATA JX10/0/

      CALL date_and_time(dtday,dttime,dtzone,idatetime)

      WRITE(6,*)
      WRITE(6,*)' Starting spectrum calculations according to SCHWINGER: '
     &  ,dttime(1:2),':',dttime(3:4),':',dttime(5:6)
      WRITE(6,*)

      GAMMA=DMYGAMMA

      IF (APERTHICK.GT.0.0D0) THEN

        APERV(2)=SIN(APERVANG)
        APERV(3)=SIN(APERHANG)
        APERV(1)=SQRT(1.0D0-(APERV(2)**2+APERV(3)**2))

      ELSE

        APERV=0.0D0
        APERHANG=0.0D0
        APERVANG=0.0D0

      ENDIF

      IF (NOBSV.GT.1) THEN
        WRITE(6,*)' '
        WRITE(6,*)
     &    '      counting from 1 to 10 for first source to show progress:'
        WRITE(6,*)' '
      ENDIF

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  '     Subroutine SPECDIP4 called to calculate spectra according to SCHWINGER'
      WRITE(LUNGFO,*)
     &  '     for found sources:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  '     photon energy cut-off parameter (SPECCUT): ',SNGL(SPECCUT)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     thickness of aperture pinhole [m]:',sngl(aperthick)
      WRITE(LUNGFO,*)'     hor. and vert. angle of aperture [rad]:',
     &  sngl(aperhang),sngl(apervang)
      WRITE(LUNGFO,*)

      JDX10=NOBSV/10
      JX10=JDX10
      IX10=1

      IF (JDX10.LT.1) JDX10=1

      DO ISOUR=1,NSOURCE

        IWARN=0
        IWARNOB=0
        DO JC=1,4
          DO IC=1,3
            SOURCEAO(IC,JC,ISOUR)=SOURCEA(IC,JC,ISOUR)
            SOURCEEO(IC,JC,ISOUR)=SOURCEE(IC,JC,ISOUR)
          ENDDO
        ENDDO
        CALL TRASOU(ISOUR)

        DO IOBSV=1,NOBSV

          DO I=1,MCO

            T=DWT(I)

            CALL WAVE_TRACK_INTER
     &        (T,XT,YT,ZT,VXT,VYT,VZT,VXP,VYP,VZP,BS,ISPLN,GAMMA)

            RX=OBSV(1,IOBSV)-XT
            RY=OBSV(2,IOBSV)-YT
            RZ=OBSV(3,IOBSV)-ZT

            IF (RX.LE.0.D0) THEN
              WRITE(LUNGFO,*)'*** ERROR IN SPECDIPA:'
              WRITE(LUNGFO,*)'BAD X-DISTANCE FROM SOURCE TO OBSERVER'
              WRITE(LUNGFO,*)'CHECK INPUT FILE'
              WRITE(6,*)'*** ERROR IN SPECDIPA:'
              WRITE(6,*)'BAD X-DISTANCE FROM SOURCE TO OBSERVER'
              WRITE(6,*)'CHECK INPUT FILE'
              STOP '--- PROGRAM WAVE ABORTED ---'
            ENDIF

            RR=RX*RX+RY*RY+RZ*RZ
            RN=SQRT(RR)
            VN=SQRT(VXT*VXT+VYT*VYT+VZT*VZT)
            VXT=VXT/VN
            VYT=VYT/VN
            VZT=VZT/VN
            OBANG=ACOS(MIN((RX*VXT+RY*VYT+RZ*VZT)/RN,1.D0))

            CALL MYBFELD(XT,YT,ZT,BX,BY,BZ,DUM,DUM,DUM)

            BS=BX*VXT+BY*VYT+BZ*VZT
            BX=BX-BS*VXT
            BY=BY-BS*VYT
            BZ=BZ-BS*VZT
            BS=SQRT(BX*BX+BY*BY+BZ*BZ)
            BRX=BY*RZ-BZ*RY
            BRY=BZ*RX-BX*RZ
            BRZ=BX*RY-BY*RX
            BRN=SQRT(BRX*BRX+BRY*BRY+BRZ*BRZ)

            IF (BS.NE.0.D0) THEN
              PSI=ACOS(MIN(BRN/BS/RN,1.D0))
     &          *SIGN(1.D0,BX*RX+BY*RY+BZ*RZ)
            ELSE
              PSI=0.D0
            ENDIF

            X0DIP(1)=XT
            Y0DIP(1)=YT
            Z0DIP(1)=ZT
            B0DIP(1)=BS

            IF (B0DIP(1).NE.0.D0) THEN
              RHODIP(1)=EMOM/B0DIP(1)/CLIGHT1
            ELSE
              WRITE(LUNGFO,*)'*** ERROR IN SPECDIP ***'
              WRITE(LUNGFO,*)'zero magnetic field for dipole ',ISOUR
              WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
              WRITE(6,*)'*** ERROR IN SPECDIPA ***'
              WRITE(6,*)'zero magnetic field for dipole ',ISOUR
              WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
              STOP
            ENDIF

            OMEGAC=1.5D0*GAMMA**3*CLIGHT1/RHODIP(1)
            ECDIP(1)=OMEGAC*HBAR1/ECHARGE1

            IF (APERTHICK.GT.0.0D0) THEN

              RXX=OBSV(1,ICBRILL)-X0DIP(ISOUR)
              RYY=OBSV(2,ICBRILL)-Y0DIP(ISOUR)
              RZZ=OBSV(3,ICBRILL)-Z0DIP(ISOUR)

              IF (RXX.LE.0.D0) THEN
                WRITE(LUNGFO,*)'*** ERROR IN SPECDIP:'
                WRITE(LUNGFO,*)'BAD X-DISTANCE FROM SOURCE TO OBSERVER'
                WRITE(LUNGFO,*)'CHECK INPUT FILE'
                WRITE(6,*)'*** ERROR IN SPECDIP:'
                WRITE(6,*)'BAD X-DISTANCE FROM SOURCE TO OBSERVER'
                WRITE(6,*)'CHECK INPUT FILE'
                STOP '--- PROGRAM WAVE ABORTED ---'
              ENDIF

              RRR=RXX*RXX+RYY*RYY+RZZ*RZZ
              RRN=SQRT(RRR)

              APERANG=ACOS((RXX*APERV(1)+RYY*APERV(2)+RZZ*APERV(3))/RRN)

              IF (IPINCIRC.EQ.0) THEN
                CALL THICKAPP(IPINCIRC,APERTHICK,PINW/2.0D0,APERANG,APERCORR)
              ELSE !IPINCIRC
                CALL THICKAPP(IPINCIRC,APERTHICK,PINR,APERANG,APERCORR)
              ENDIF

            ELSE !APERTHICK

              APERCORR=1.0D0

            ENDIF !APERTHICK

            DO IFREQ=1,NFREQ

              Y=FREQ(IFREQ)/ECDIP(1)

              ILIOBFR=ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1))
              IOBFR=IOBSV+NOBSV*(IFREQ-1)

              IF (XT.GE.XIANF.AND.XT.LE.XIEND) THEN
                speck=
     &            DFDTDP(Y,PSI,GAMMA,DMYCUR,BANWID,PAR,PER,POWR)/RR
     &            *APERCORR
                IF (SPECCUT.NE.0.D0.AND.Y.GT.SPECCUT) THEN
                  speck=0.0d0
                  PAR=0.D0
                  PER=0.D0
                  powr=0.0d0
                ENDIF   !SPECCUT

                SPEC(ILIOBFR)=SPEC(ILIOBFR)+speck

                speck=POWR/RHODIP(1)/RR*APERCORR
                SPECPOW(ISOUR+NSOURCE*(IOBSV-1))=
     &            SPECPOW(ISOUR+NSOURCE*(IOBSV-1))+speck

                EXPOM=(1.D0,0.D0)

                IF (IPIN.NE.0) THEN
                  DX2=RX*RX
                  DZY2=RZ*RZ+RY*RY
                  OM=FREQ(IFREQ)/(HBAREV1*CLIGHT1)

C       TO MAKE SURE THAT TAYLOR-EXPANSION IS VALID

                  IF (DZY2.GT.0.01D0*DX2) THEN
                    WRITE(LUNGFO,*)
     &                '*** ERROR IN SPECDIPA: OBSERVATION ANGLE TO LARGE ***'
                    WRITE(LUNGFO,*)'DECREASE SIZE OF PINHOLE OR WGWINFC ...'
                    WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
                    WRITE(6,*)
     &                '*** ERROR IN SPECDIPA: OBSERVATION ANGLE TO LARGE ***'
                    WRITE(6,*)'DECREASE SIZE OF PINHOLE OR WGWINFC ...'
                    WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
                    STOP
                  ENDIF   !(DZY2.GT.0.01D0*DX2)

                  EPS(1)=DZY2/DX2
                  DO IEPS=2,6
                    EPS(IEPS)=EPS(IEPS-1)*EPS(1)
                  ENDDO !IEPS

c        TAYLOR-EXPANSION DONE WITH REDUCE
c       IN "WTAY1.RED";
c       on rounded;
c       on numval;
c       precision 13;
c       F:=SQRT(1+EPS);
c       DR:=TAY1(F,EPS,6);
c       ON FORT;
c       OUT "RED.FOR";
c       DR;
c       SHUT "RED.FOR";
C ans is actually reduce by 1.0 to avoid large overall phase

                  ans=-0.0205078125D0*eps(6)+0.02734375D0*eps(5)
     &              -0.0390625D0*eps(4)+
     &              0.0625D0*eps(3)-0.125D0*eps(2)+0.5D0*eps(1)

                  DRRED=-DABS(RX*ANS)
                  EXPOM=CDEXP(DCMPLX(0.D0,DRRED*OM))

                ENDIF !IPIN

                AFREQ(1)=(0.D0,0.D0)
                AFREQ(2)=afreq(2)+DCMPLX(0.D0,-SQRT(PER/RR/SPECNOR))*EXPOM
     &            *APERCORR
                AFREQ(3)=afreq(3)+DCMPLX(SQRT(PAR/RR/SPECNOR),0.D0)*EXPOM*SIGN(1.D0,PSI)
     &            *APERCORR

                IF (IWARNOB.EQ.0.AND.OBANGMN.GT.WGWINFC/GAMMA) THEN
                  WRITE(LUNGFO,*)'*** WARNING IN SPECDIPA:'
                  WRITE(LUNGFO,*)'problems finding tangent point'
                  WRITE(LUNGFO,*)'source number ',ISOUR
                  WRITE(LUNGFO,*)
     &              'maybe low WBL0CUT or WGWINFC in namelist COLLIN causes problems'
                  WRITE(LUNGFO,*)
     &              '... or turn option IWIGGLER off and tune parameters'
                  WRITE(LUNGFO,*)
     &              'ISPECDIP, WBL0CUT ... by hand'
                  WRITE(LUNGFO,*)
                  WRITE(LUNGFO,*)' *** Photon flux set to zero! ***'
                  WRITE(LUNGFO,*)
                  WRITE(6,*)'*** WARNING IN SPECDIPA:'
                  WRITE(6,*)'problems finding tangent point'
                  WRITE(6,*)'source number ',ISOUR
                  WRITE(6,*)
     &              'maybe low WBL0CUT or WGWINFC in namelist COLLIN causes problems'
                  WRITE(6,*)
     &              '... or turn option IWIGGLER off and tune parameters'
                  WRITE(6,*)
     &              'ISPECDIP, WBL0CUT ... by hand'
                  WRITE(6,*)
                  WRITE(6,*)' *** Photon flux set to zero! ***'
                  WRITE(6,*)
                  IWARNOB=1
                ENDIF

                REAIMA(1,1,IOBFR)=DREAL(AFREQ(1))
                REAIMA(1,2,IOBFR)=DIMAG(AFREQ(1))
                REAIMA(2,1,IOBFR)=DREAL(AFREQ(2))
                REAIMA(2,2,IOBFR)=DIMAG(AFREQ(2))
                REAIMA(3,1,IOBFR)=DREAL(AFREQ(3))
                REAIMA(3,2,IOBFR)=DIMAG(AFREQ(3))
                IF (IPOLA.NE.0) THEN
                  APOL=
     &              AFREQ(1)*CONJG(VPOLA(1))
     &              +AFREQ(2)*CONJG(VPOLA(2))
     &              +AFREQ(3)*CONJG(VPOLA(3))
                  SPEC(ILIOBFR)=
     &              DREAL(APOL*CONJG(APOL))*SPECNOR
                ENDIF !IPOLA

                IF (ISTOKES.NE.0) THEN

                  APOLH=
     &              AFREQ(1)*CONJG(VSTOKES(1,1))
     &              +AFREQ(2)*CONJG(VSTOKES(1,2))
     &              +AFREQ(3)*CONJG(VSTOKES(1,3))

                  APOLR=
     &              AFREQ(1)*CONJG(VSTOKES(2,1))
     &              +AFREQ(2)*CONJG(VSTOKES(2,2))
     &              +AFREQ(3)*CONJG(VSTOKES(2,3))

                  APOLL=
     &              AFREQ(1)*CONJG(VSTOKES(3,1))
     &              +AFREQ(2)*CONJG(VSTOKES(3,2))
     &              +AFREQ(3)*CONJG(VSTOKES(3,3))

                  APOL45=
     &              AFREQ(1)*CONJG(VSTOKES(4,1))
     &              +AFREQ(2)*CONJG(VSTOKES(4,2))
     &              +AFREQ(3)*CONJG(VSTOKES(4,3))

                  STOK1=
     &              REAL(APOLR*CONJG(APOLR))+
     &              REAL(APOLL*CONJG(APOLL))

                  STOK2=-STOK1+
     &              2.*REAL(APOLH*CONJG(APOLH))

                  STOK3=
     &              2.*REAL(APOL45*CONJG(APOL45))-
     &              STOK1

                  STOK4=
     &              REAL(APOLR*CONJG(APOLR))-
     &              REAL(APOLL*CONJG(APOLL))

                  STOKES(1,IOBFR)=STOKES(1,IOBFR)+
     &              STOK1*SPECNOR

                  STOKES(2,IOBFR)=STOKES(2,IOBFR)+
     &              STOK2*SPECNOR

                  STOKES(3,IOBFR)=STOKES(3,IOBFR)+
     &              STOK3*SPECNOR

                  STOKES(4,IOBFR)=STOKES(4,IOBFR)+
     &              STOK4*SPECNOR

                ENDIF !ISTOKES

              ENDIF   !XIANF

              SPECTOT(IOBFR)=SPECTOT(IOBFR)
     &          +SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))

            ENDDO   !IFREQ

            IF (NOBSV.GT.1.AND.ISOUR.EQ.1.AND.IOBSV.EQ.JX10.AND.IX10.LE.10) THEN
              JX10=JX10+JDX10
              CALL date_and_time(dtday,dttime,dtzone,idatetime)
              WRITE(6,*)' ',IX10,' ',dttime(1:2),':',dttime(3:4),':',dttime(5:6)
              IX10=IX10+1
            ENDIF

          enddo !mco

        ENDDO   !IOBSV

        IF (ISPECDIP.LE.-3.and.isour.eq.nsource) THEN
          DEALLOCATE(DWT)
          DEALLOCATE(DWX)
          DEALLOCATE(DWX2P)
          DEALLOCATE(DWB)
          DEALLOCATE(DWB2P)
          DEALLOCATE(DWY)
          DEALLOCATE(DWY2P)
          DEALLOCATE(DWZ)
          DEALLOCATE(DWZ2P)
          DEALLOCATE(tragam)
        ENDIF

        IF (ISOUR.EQ.1.and.nsource.gt.1) THEN
          WRITE(6,*)' '
          WRITE(6,*)' '
          WRITE(6,*)' sources treated so far:'
          WRITE(6,*)' '
        ENDIF

        CALL date_and_time(dtday,dttime,dtzone,idatetime)
        WRITE(6,2000)ISOUR,NSOURCE,dttime(1:2),dttime(3:4),dttime(5:6)
2000    FORMAT(10X,I4,' of',I4,2X,A,':',A,':',A)

      ENDDO   !ISOUR

      WRITE(6,*)' '

      RETURN
      END
