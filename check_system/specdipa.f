*CMZ :  4.00/17 04/10/2022  08.10.22  by  Michael Scheer
*CMZ :  4.00/15 13/02/2022  16.33.42  by  Michael Scheer
*CMZ :  4.00/11 15/06/2021  10.34.53  by  Michael Scheer
*CMZ :  3.08/01 02/04/2019  15.33.15  by  Michael Scheer
*CMZ :  3.02/03 07/11/2014  15.53.20  by  Michael Scheer
*CMZ :  3.02/00 18/09/2014  14.19.02  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  14.04.56  by  Michael Scheer
*CMZ :  2.69/00 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.66/04 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.66/03 12/11/2009  16.27.11  by  Michael Scheer
*CMZ :  2.63/05 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.57/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.52/03 08/07/2004  13.35.44  by  Michael Scheer
*CMZ :  2.44/01 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.41/09 14/08/2002  17.16.34  by  Michael Scheer
*CMZ :  2.41/04 21/03/2002  12.47.15  by  Michael Scheer
*CMZ :  2.41/02 21/03/2002  12.39.45  by  Michael Scheer
*CMZ :  2.41/01 20/03/2002  19.34.03  by  Michael Scheer
*CMZ :  2.41/00 20/03/2002  19.23.14  by  Michael Scheer
*CMZ :  2.34/01 01/06/2001  15.33.00  by  Michael Scheer
*CMZ :  2.33/02 03/05/2001  17.18.32  by  Michael Scheer
*CMZ :  2.33/01 03/05/2001  12.21.08  by  Michael Scheer
*CMZ :  2.33/00 03/05/2001  10.55.42  by  Michael Scheer
*CMZ :  2.32/04 26/04/2001  12.24.00  by  Michael Scheer
*CMZ :  2.32/01 25/04/2001  19.08.58  by  Michael Scheer
*CMZ :  2.31/01 25/04/2001  17.10.34  by  Michael Scheer
*CMZ :  2.31/00 24/04/2001  14.48.24  by  Michael Scheer
*CMZ :  2.20/01 12/12/2000  10.01.31  by  Michael Scheer
*CMZ :  2.16/08 31/10/2000  17.29.14  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.13/03 10/01/2000  17.17.14  by  Michael Scheer
*CMZ :  2.00/00 15/12/98  14.29.27  by  Michael Scheer
*CMZ : 00.02/02 15/01/97  15.18.41  by  Michael Scheer
*CMZ : 00.00/07 18/05/94  14.54.15  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.54.12  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.24  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE SPECDIPA
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
*KEEP,debugwave.
      include 'debugwave.cmn'
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
     &  GAMMA

      REAL*8 STOK1,STOK2,STOK3,STOK4

      COMPLEX*16 APOL,AFREQ(3),EXPOM
      COMPLEX*16 APOLH,APOLR,APOLL,APOL45

      DATA JX10/0/

      IF (ISPECDIP.EQ.-4) THEN
        call specdip4
        return
      endif

      CALL date_and_time(dtday,dttime,dtzone,idatetime)

      if (ispecdip.ne.0) then
        WRITE(6,*)
        WRITE(6,*)' Starting spectrum calculations according to SCHWINGER: '
     &    ,dttime(1:2),':',dttime(3:4),':',dttime(5:6)
        WRITE(6,*)
      endif

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

      if (ispecdip.ne.0) then
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     Subroutine SPECDIPA called to calculate spectra according to SCHWINGER'
        WRITE(LUNGFO,*)
     &    '     for found sources:'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     photon energy cut-off parameter (SPECCUT): ',SNGL(SPECCUT)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     thickness of aperture pinhole [m]:',sngl(aperthick)
        WRITE(LUNGFO,*)'     hor. and vert. angle of aperture [rad]:',
     &    sngl(aperhang),sngl(apervang)
        WRITE(LUNGFO,*)
      endif

      JDX10=NOBSV/10
      JX10=JDX10
      IX10=1

      IF (JDX10.LT.1) JDX10=1

      IF (ISPECDIP.EQ.-1) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '      center of found source is taken for as dipole source'
        WRITE(LUNGFO,*)
      ELSE IF (ISPECDIP.EQ.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '      Points of primary trajectory with lowest radiation angle'
        WRITE(LUNGFO,*)
     &    '      to observation points are taken to estimate location of source points.'
        WRITE(LUNGFO,*)
      ELSE IF (ISPECDIP.EQ.-2) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '      Points of primary trajectory with lowest radiation angle'
        WRITE(LUNGFO,*)
     &    '      to observation points are taken as dipole source.'
        WRITE(LUNGFO,*)
      ELSE
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '      Source from primary trajectory is retracked with high resolution.'
        WRITE(LUNGFO,*)
     &    '      Points with lowest radiation angle to observation points'
        WRITE(LUNGFO,*)
     &    '      are taken as dipole sources.'
        WRITE(LUNGFO,*)
      ENDIF

      WRITE(LUNGFO,*)
     &  '     number of source, number of observation point,'
      WRITE(LUNGFO,*)
     &  '     field, bending radius ,Ec [eV], cut-off (due to SPECCUT),'
      WRITE(LUNGFO,*)
     &  '      and X,Y,Z-Pos. [m] for selected observation point:'
      WRITE(LUNGFO,*)

      DO ISOUR=1,NSOURCE

        IWARN=0
        IWARNOB=0

        IF (ISPECDIP.LE.-3) THEN
          DO JC=1,4
            DO IC=1,3
              SOURCEAO(IC,JC,ISOUR)=SOURCEA(IC,JC,ISOUR)
              SOURCEEO(IC,JC,ISOUR)=SOURCEE(IC,JC,ISOUR)
            ENDDO
          ENDDO
          CALL TRASOU(ISOUR)
        ENDIF  !ISPECDIP

        DO IOBSV=1,NOBSV

          IF (ISPECDIP.EQ.-2.or.ispecdip.eq.0) THEN

            OBANGMN=1.D30
            IF (ISOURAE(2,ISOUR)-ISOURAE(1,ISOUR).LT.3) THEN
              ITANG=ISOURAE(1,ISOUR)+1
            ELSE  !ISOURAE
              DO I=ISOURAE(1,ISOUR)+1,ISOURAE(2,ISOUR)-1
                XT=WTRA(1,1,I)
                YT=WTRA(2,1,I)
                ZT=WTRA(3,1,I)
                VXT=WTRA(1,2,I)
                VYT=WTRA(2,2,I)
                VZT=WTRA(3,2,I)
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
                OBANG=ACOS(MIN((RX*VXT+RY*VYT+RZ*VZT)/RN/VN,1.D0))
                IF (ABS(OBANG).LT.OBANGMN) THEN
                  ITANG=I
                  OBANGMN=OBANG
                ENDIF   !(ABS(OBANG).LT.OBANGMN)
              ENDDO   !ISOURAE
            ENDIF   !ISOURAE

            IF (ITANG.LT.2) ITANG=2
            IF (ITANG.GE.NCO) ITANG=NCO-1
            DO J=1,3
              I=ITANG+J-2
              XT=WTRA(1,1,I)
              YT=WTRA(2,1,I)
              ZT=WTRA(3,1,I)
              VXT=WTRA(1,2,I)
              VYT=WTRA(2,2,I)
              VZT=WTRA(3,2,I)
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
              OBANG=ACOS(MIN((RX*VXT+RY*VYT+RZ*VZT)/RN/VN,1.D0))
              ANG3(J)=OBANG
              XT3(J)=XT
              YT3(J)=YT
              ZT3(J)=ZT
              VXT3(J)=VXT
              VYT3(J)=VYT
              VZT3(J)=VZT
            ENDDO   !J=1,3

            CALL UTIL_PARABEL(XT3,ANG3,A,DUM3,XT,OBANG,IFAIL)

            IF (XT.GT.XT3(3).OR.XT.LT.XT3(1)) THEN

              IF (IWARN.EQ.0) THEN
                WRITE(LUNGFO,*)'*** WARNING IN SPECDIPA:'
                WRITE(LUNGFO,*)'problems finding tangent point'
                WRITE(LUNGFO,*)'source number ',ISOUR
                WRITE(LUNGFO,*)
     &            'maybe low WBL0CUT in namelist COLLIN causes problems'
                WRITE(LUNGFO,*)'presumably you may ignore this warning'
                WRITE(LUNGFO,*)
     &            '... or turn option IWIGGLER off and tune parameters'
                WRITE(LUNGFO,*)
     &            'ISPECDIP, WBL0CUT ... by hand'
                WRITE(6,*)'*** WARNING IN SPECDIPA:'
                WRITE(6,*)'problems finding tangent point'
                WRITE(6,*)'source number ',ISOUR
                WRITE(6,*)
     &            'maybe low WBL0CUT in namelist COLLIN causes problems'
                WRITE(6,*)'presumably you may ignore this warning'
                WRITE(6,*)
     &            '... or turn option IWIGGLER off and tune parameters'
                WRITE(6,*)
     &            'ISPECDIP, WBL0CUT ... by hand'
                IWARN=1
              ENDIF   !IWARN

              XT=XT3(2)

            ENDIF   !XT3

            IF (IFAIL.NE.0) THEN
              WRITE(LUNGFO,*)'*** WARNING IN SPECDIPA:'
              WRITE(LUNGFO,*)'BAD RETURN OF UTIL_PARABEL'
              WRITE(LUNGFO,*)'TRY DIFFERENT ISPECDIP OR OTHER SETTING'
              WRITE(6,*)'*** WARNING IN SPECDIPA:'
              WRITE(6,*)'BAD RETURN OF UTIL_PARABEL'
              WRITE(6,*)'TRY DIFFERENT ISPECDIP OR OTHER SETTING'
            ENDIF   !IFAIL

            CALL UTIL_PARABEL(XT3,YT3,A,DUM3,DUM,DUM,IFAIL)

            IF (IFAIL.NE.0) THEN
              WRITE(LUNGFO,*)'*** WARNING IN SPECDIPA:'
              WRITE(LUNGFO,*)'BAD RETURN OF UTIL_PARABEL'
              WRITE(LUNGFO,*)'TRY DIFFERENT ISPECDIP OR OTHER SETTING'
              WRITE(6,*)'*** WARNING IN SPECDIPA:'
              WRITE(6,*)'BAD RETURN OF UTIL_PARABEL'
              WRITE(6,*)'TRY DIFFERENT ISPECDIP OR OTHER SETTING'
            ENDIF   !IFAIL

            YT=A(1)+A(2)*XT+A(3)*XT*XT

            CALL UTIL_PARABEL(XT3,ZT3,A,DUM3,DUM,DUM,IFAIL)

            IF (IFAIL.NE.0) THEN
              WRITE(LUNGFO,*)'*** WARNING IN SPECDIPA:'
              WRITE(LUNGFO,*)'BAD RETURN OF UTIL_PARABEL'
              WRITE(LUNGFO,*)'TRY DIFFERENT ISPECDIP OR OTHER SETTING'
              WRITE(6,*)'*** WARNING IN SPECDIPA:'
              WRITE(6,*)'BAD RETURN OF UTIL_PARABEL'
              WRITE(6,*)'TRY DIFFERENT ISPECDIP OR OTHER SETTING'
            ENDIF   !IFAIL

            ZT=A(1)+A(2)*XT+A(3)*XT*XT

            CALL UTIL_PARABEL(XT3,VXT3,A,DUM3,DUM,DUM,IFAIL)

            IF (IFAIL.NE.0) THEN
              WRITE(LUNGFO,*)'*** WARNING IN SPECDIPA:'
              WRITE(LUNGFO,*)'BAD RETURN OF UTIL_PARABEL'
              WRITE(LUNGFO,*)'TRY DIFFERENT ISPECDIP OR OTHER SETTING'
              WRITE(6,*)'*** WARNING IN SPECDIPA:'
              WRITE(6,*)'BAD RETURN OF UTIL_PARABEL'
              WRITE(6,*)'TRY DIFFERENT ISPECDIP OR OTHER SETTING'
            ENDIF   !IFAIL

            VXT=A(1)+A(2)*XT+A(3)*XT*XT

            CALL UTIL_PARABEL(XT3,VYT3,A,DUM3,DUM,DUM,IFAIL)

            IF (IFAIL.NE.0) THEN
              WRITE(LUNGFO,*)'*** WARNING IN SPECDIPA:'
              WRITE(LUNGFO,*)'BAD RETURN OF UTIL_PARABEL'
              WRITE(LUNGFO,*)'TRY DIFFERENT ISPECDIP OR OTHER SETTING'
              WRITE(6,*)'*** WARNING IN SPECDIPA:'
              WRITE(6,*)'BAD RETURN OF UTIL_PARABEL'
              WRITE(6,*)'TRY DIFFERENT ISPECDIP OR OTHER SETTING'
            ENDIF   !IFAIL

            VYT=A(1)+A(2)*XT+A(3)*XT*XT

            CALL UTIL_PARABEL(XT3,VZT3,A,DUM3,DUM,DUM,IFAIL)

            IF (IFAIL.NE.0) THEN
              WRITE(LUNGFO,*)'*** WARNING IN SPECDIPA:'
              WRITE(LUNGFO,*)'BAD RETURN OF UTIL_PARABEL'
              WRITE(LUNGFO,*)'TRY DIFFERENT ISPECDIP OR OTHER SETTING'
              WRITE(6,*)'*** WARNING IN SPECDIPA:'
              WRITE(6,*)'BAD RETURN OF UTIL_PARABEL'
              WRITE(6,*)'TRY DIFFERENT ISPECDIP OR OTHER SETTING'
            ENDIF   !IFAIL

            VZT=A(1)+A(2)*XT+A(3)*XT*XT

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

            CALL MYBFELD(XT,YT,ZT,BX,BY,BZ,DUM,DUM,DUM)

            VN=SQRT(VXT*VXT+VYT*VYT+VZT*VZT)
            VXT=VXT/VN
            VYT=VYT/VN
            VZT=VZT/VN
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

          ELSE IF (ISPECDIP.LE.-3) THEN

            OBANGMN=1.D30
            ISPLN=0

            DO I=1,MCO

              T=DWT(I)

              CALL WAVE_TRACK_INTER
     &          (T,XT,YT,ZT,VXT,VYT,VZT,VXP,VYP,VZP,BS,ISPLN,GAMMA)

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

              IF (OBANG.LT.OBANGMN) THEN

                OBANGMN=OBANG

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
     &              *SIGN(1.D0,BX*RX+BY*RY+BZ*RZ)
                ELSE
                  PSI=0.D0
                ENDIF

                X0DIP(1)=XT
                Y0DIP(1)=YT
                Z0DIP(1)=ZT
                B0DIP(1)=BS
                RR3=RR

              ENDIF  !(ABS(OBANG).LT.OBANGMN)

            ENDDO   !MCO

            RR=RR3

          ELSE IF (ISPECDIP.EQ.-1) THEN

            XT=SOURCEN(1,1,ISOUR)
            YT=SOURCEN(2,1,ISOUR)
            ZT=SOURCEN(3,1,ISOUR)

            VXT=SOURCEN(1,2,ISOUR)
            VYT=SOURCEN(2,2,ISOUR)
            VZT=SOURCEN(3,2,ISOUR)
            VN=SQRT(VXT*VXT+VYT*VYT+VZT*VZT)
            VXT=VXT/VN
            VYT=VYT/VN
            VZT=VZT/VN

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

          ENDIF   !ISPECDIP

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

c          IF (IOBSV.EQ.ICBRILL) THEN
          WRITE (LUNGFO,*)'      ',isour,iobsv
          WRITE (LUNGFO,*)'      ',
     &      SNGL(B0DIP(1)),SNGL(RHODIP(1)),SNGL(ECDIP(1))
     &      ,SNGL(ECDIP(1)*SPECCUT)
          WRITE(LUNGFO,*)'                   ',
     &      SNGL(X0DIP(1)),SNGL(Y0DIP(1)),SNGL(Z0DIP(1))
c          ENDIF !ICBRILL

          schwingercen(1,iobsv,isour)=x0dip(1)
          schwingercen(2,iobsv,isour)=y0dip(1)
          schwingercen(3,iobsv,isour)=z0dip(1)
          schwingercen(4,iobsv,isour)=b0dip(1)

          if (ispecdip.eq.0) cycle

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
              SPEC(ILIOBFR)=
     &          DFDTDP(Y,PSI,GAMMA,DMYCUR,BANWID,PAR,PER,POWR)/RR
     &          *APERCORR
              IF (SPECCUT.NE.0.D0.AND.Y.GT.SPECCUT) THEN
                SPEC(ILIOBFR)=0.D0
                PAR=0.D0
                PER=0.D0
                powr=0.0d0
              ENDIF   !SPECCUT

              SPECPOW(ISOUR+NSOURCE*(IOBSV-1))=POWR/RHODIP(1)/RR
     &          *APERCORR

              EXPOM=(1.D0,0.D0)

              IF (IPIN.NE.0) THEN
                DX2=RX*RX
                DZY2=RZ*RZ+RY*RY
                OM=FREQ(IFREQ)/(HBAREV1*CLIGHT1)

C       TO MAKE SURE THAT TAYLOR-EXPANSION IS VALID

                IF (DZY2.GT.0.01D0*DX2) THEN
                  WRITE(LUNGFO,*)
     &              '*** ERROR IN SPECDIPA: OBSERVATION ANGLE TO LARGE ***'
                  WRITE(LUNGFO,*)'DECREASE SIZE OF PINHOLE OR WGWINFC ...'
                  WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
                  WRITE(6,*)
     &              '*** ERROR IN SPECDIPA: OBSERVATION ANGLE TO LARGE ***'
                  WRITE(6,*)'DECREASE SIZE OF PINHOLE OR WGWINFC ...'
                  WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
                  STOP
                ENDIF     !(DZY2.GT.0.01D0*DX2)

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
     &            -0.0390625D0*eps(4)+
     &            0.0625D0*eps(3)-0.125D0*eps(2)+0.5D0*eps(1)

                DRRED=-DABS(RX*ANS)
                EXPOM=CDEXP(DCMPLX(0.D0,DRRED*OM))

              ENDIF !IPIN

              AFREQ(1)=(0.D0,0.D0)
              AFREQ(2)=DCMPLX(0.D0,-SQRT(PER/RR/SPECNOR))*EXPOM
     &          *APERCORR
              AFREQ(3)=DCMPLX(SQRT(PAR/RR/SPECNOR),0.D0)*EXPOM*SIGN(1.D0,PSI)
     &          *APERCORR

              IF (IWARNOB.EQ.0.AND.OBANGMN.GT.WGWINFC/GAMMA) THEN
                WRITE(LUNGFO,*)'*** WARNING IN SPECDIPA:'
                WRITE(LUNGFO,*)'problems finding tangent point'
                WRITE(LUNGFO,*)'source number ',ISOUR
                WRITE(LUNGFO,*)
     &            'maybe low WBL0CUT or WGWINFC in namelist COLLIN causes problems'
                WRITE(LUNGFO,*)
     &            '... or turn option IWIGGLER off and tune parameters'
                WRITE(LUNGFO,*)
     &            'ISPECDIP, WBL0CUT ... by hand'
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)' *** Photon flux set to zero! ***'
                WRITE(LUNGFO,*)
                WRITE(6,*)'*** WARNING IN SPECDIPA:'
                WRITE(6,*)'problems finding tangent point'
                WRITE(6,*)'source number ',ISOUR
                WRITE(6,*)
     &            'maybe low WBL0CUT or WGWINFC in namelist COLLIN causes problems'
                WRITE(6,*)
     &            '... or turn option IWIGGLER off and tune parameters'
                WRITE(6,*)
     &            'ISPECDIP, WBL0CUT ... by hand'
                WRITE(6,*)
                WRITE(6,*)' *** Photon flux set to zero! ***'
                WRITE(6,*)
                IWARNOB=1
              ENDIF

              IF (OBANGMN.GT.WGWINFC/GAMMA) THEN
                AFREQ(1)=(0.D0,0.D0)
                AFREQ(2)=(0.D0,0.D0)
                AFREQ(3)=(0.D0,0.D0)
              ENDIF

              REAIMA(1,1,IOBFR)=DREAL(AFREQ(1))
              REAIMA(1,2,IOBFR)=DIMAG(AFREQ(1))
              REAIMA(2,1,IOBFR)=DREAL(AFREQ(2))
              REAIMA(2,2,IOBFR)=DIMAG(AFREQ(2))
              REAIMA(3,1,IOBFR)=DREAL(AFREQ(3))
              REAIMA(3,2,IOBFR)=DIMAG(AFREQ(3))
              IF (IPOLA.NE.0) THEN
                APOL=
     &            AFREQ(1)*CONJG(VPOLA(1))
     &            +AFREQ(2)*CONJG(VPOLA(2))
     &            +AFREQ(3)*CONJG(VPOLA(3))
                SPEC(ILIOBFR)=
     &            DREAL(APOL*CONJG(APOL))*SPECNOR
              ENDIF !IPOLA

              IF (ISTOKES.NE.0) THEN

                APOLH=
     &            AFREQ(1)*CONJG(VSTOKES(1,1))
     &            +AFREQ(2)*CONJG(VSTOKES(1,2))
     &            +AFREQ(3)*CONJG(VSTOKES(1,3))

                APOLR=
     &            AFREQ(1)*CONJG(VSTOKES(2,1))
     &            +AFREQ(2)*CONJG(VSTOKES(2,2))
     &            +AFREQ(3)*CONJG(VSTOKES(2,3))

                APOLL=
     &            AFREQ(1)*CONJG(VSTOKES(3,1))
     &            +AFREQ(2)*CONJG(VSTOKES(3,2))
     &            +AFREQ(3)*CONJG(VSTOKES(3,3))

                APOL45=
     &            AFREQ(1)*CONJG(VSTOKES(4,1))
     &            +AFREQ(2)*CONJG(VSTOKES(4,2))
     &            +AFREQ(3)*CONJG(VSTOKES(4,3))

                STOK1=
     &            APOLR*CONJG(APOLR)+
     &            APOLL*CONJG(APOLL)

                STOK2=-STOK1+
     &            2.0d0*APOLH*CONJG(APOLH)

                STOK3=
     &            2.0d0*APOL45*CONJG(APOL45)-
     &            STOK1

                STOK4=
     &            APOLR*CONJG(APOLR)-
     &            APOLL*CONJG(APOLL)

                if (abs(stok1)*specnor.gt.1.0D-30)
     &            STOKES(1,IOBFR)=STOKES(1,IOBFR)+
     &            STOK1*SPECNOR

                if (abs(stok2)*specnor.gt.1.0D-30)
     &            STOKES(2,IOBFR)=STOKES(2,IOBFR)+
     &            STOK2*SPECNOR

                if (abs(stok3)*specnor.gt.1.0D-30)
     &            STOKES(3,IOBFR)=STOKES(3,IOBFR)+
     &            STOK3*SPECNOR

                if (abs(stok4)*specnor.gt.1.0D-30)
     &           STOKES(4,IOBFR)=STOKES(4,IOBFR)+
     &            STOK4*SPECNOR

              ENDIF !ISTOKES

            ENDIF   !XIANF

            SPECTOT(IOBFR)=SPECTOT(IOBFR)
     &        +SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))

          ENDDO   !IFREQ

          IF (NOBSV.GT.1.AND.ISOUR.EQ.1.AND.IOBSV.EQ.JX10.AND.IX10.LE.10) THEN
            JX10=JX10+JDX10
            CALL date_and_time(dtday,dttime,dtzone,idatetime)
            WRITE(6,*)' ',IX10,' ',dttime(1:2),':',dttime(3:4),':',dttime(5:6)
            IX10=IX10+1
          ENDIF

        ENDDO   !IOBSV

        if (ispecdip.eq.0) cycle

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
