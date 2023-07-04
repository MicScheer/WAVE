*CMZ :  2.63/03 02/06/2009  16.19.23  by  Michael Scheer
*CMZ :  2.61/04 29/03/2007  16.12.06  by  Michael Scheer
*CMZ :  2.57/05 14/12/2006  10.19.14  by  Michael Scheer
*CMZ :  2.52/09 29/10/2004  11.46.00  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  1.03/06 06/08/98  18.01.12  by  Michael Scheer
*CMZ :  1.02/00 18/12/97  11.55.12  by  Michael Scheer
*CMZ :  1.01/00 28/10/97  18.46.16  by  Michael Scheer
*CMZ :  1.00/00 06/06/97  18.18.12  by  Michael Scheer
*-- Author :    Michael Scheer   05/06/97

      SUBROUTINE MODULATOR(GRARAD)
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

C--- INPUT REC-MODULATOR VIA THIS ROUTINE AT THE END OF SR REC_INIT

      IMPLICIT NONE

      EXTERNAL DCOSD,DSIND
      DOUBLE PRECISION DCOSD,DSIND

*KEEP,klotz.
      include 'klotz.cmn'
*KEEP,modulator.
      include 'modulator.cmn'
*KEND.

      DOUBLE PRECISION XMOD(NMAGMODP*NSLICEP)
     &  ,YMOD(NMAGMODP*NSLICEP)
     &  ,ZMOD(NMAGMODP*NSLICEP)
     &  ,DXMOD(NMAGMODP*NSLICEP)
     &  ,DYMOD(NMAGMODP*NSLICEP)
     &  ,DZMOD(NMAGMODP*NSLICEP)
     &  ,THEMOD(NMAGMODP*NSLICEP)
     &  ,PHIMOD(NMAGMODP*NSLICEP)
     &  ,BCBCMOD(NMAGMODP*NSLICEP)

      DOUBLE PRECISION Y,DLY,Y2,Y1,BCCOS,BCSIN,RHO2,GRARAD

      INTEGER IMAG,I,KMAG,LMAG,ICAL,IREAD,JMAG

      DATA ICAL/0/

      IF (NMAGMOD.LT.0) THEN
        IREAD=1
        OPEN(UNIT=99,FILE='therot.dat',STATUS='OLD')
        read(99,*)nmagmod
      ELSE
        IREAD=0
      ENDIF

      IF (NSLICE.GT.NSLICEP) THEN
        WRITE(6,*) '*** DIMENSION NSLICEP EXCEEDED ***'
        STOP
      ENDIF

      IF (NMAGMOD.GT.NMAGMODP) THEN
        WRITE(6,*)'*** DIMENSION NMAGMODP EXCEEDED ***'
        STOP
      ENDIF

      IF (2*NMAGMOD*NSLICE+IMAGTOT.GT.NKLOTZ) THEN
        WRITE(6,*)'*** ERROR IN MODULATOR: DIMENSION NKLOTZ EXCEEDED ***'
        STOP
      ENDIF

      IF (ICAL.EQ.0) THEN
        IF (IREAD.GT.0) THEN
          DO IMAG=1,NMAGMOD
            READ(99,*)RADIMOD(IMAG),ZLENMOD(IMAG),
     &        CENMODX(IMAG),CENMODY(IMAG),CENMODZ(IMAG),
     &        THEROT(IMAG),BCMOD(IMAG)
          ENDDO
          CLOSE(99)
        ENDIF !IREAD

        IF (ITHEMSYM.EQ.1.OR.ITHEMSYM.EQ.2) THEN

          IF (NMAGMOD*2*ITHEMSYM.GT.NMAGMODP) THEN
            WRITE(6,*)'*** DIMENSION NMAGMODP EXCEEDED ***'
            STOP
          ENDIF

          IF (2*NMAGMOD*2*ITHEMSYM*NSLICE+IMAGTOT.GT.NKLOTZ) THEN
            WRITE(6,*)'*** ERROR IN MODULATOR: DIMENSION NKLOTZ EXCEEDED ***'
            STOP
          ENDIF

          IF (ITHEMSYM.EQ.1) THEN !cerror? 2jun09

            DO IMAG=1,NMAGMOD
              JMAG=IMAG+NMAGMOD
              RADIMOD(JMAG)=RADIMOD(IMAG)
              ZLENMOD(JMAG)=ZLENMOD(IMAG)
              CENMODX(JMAG)=CENMODX(IMAG)
              CENMODY(JMAG)=-CENMODY(IMAG)
              CENMODZ(JMAG)=CENMODZ(IMAG)
              THEROT(JMAG)=THEROT(IMAG)
              BCMOD(JMAG)=BCMOD(IMAG)
            ENDDO !IMAG=1,NMAGMOD

          NMAGMOD=NMAGMOD*2

          ENDIF

          IF (ITHEMSYM.EQ.2) THEN

            DO IMAG=1,NMAGMOD
              JMAG=IMAG+NMAGMOD
              RADIMOD(JMAG)=RADIMOD(IMAG)
              ZLENMOD(JMAG)=ZLENMOD(IMAG)
              CENMODX(JMAG)=-CENMODX(IMAG)
              CENMODY(JMAG)=CENMODY(IMAG)
              CENMODZ(JMAG)=CENMODZ(IMAG)
              THEROT(JMAG)=THEROT(IMAG)
              BCMOD(JMAG)=BCMOD(IMAG)
            ENDDO !IMAG=1,NMAGMOD

            NMAGMOD=NMAGMOD*2

          ENDIF !ITHESMSYM

        ENDIF !ITHESMSYM

        IF (ITHESYML.EQ.-1) THEN
          DO IMAG=1,NMAGMOD
            IF (CENMODY(IMAG).LT.0.0D0) THEROT(IMAG)=-THEROT(IMAG)
          ENDDO   !IMAG=1,NMAGMOD
        ENDIF !ITHESYMU

        IF (ITHESYMD.EQ.-1) THEN
          DO IMAG=1,NMAGMOD
            IF (CENMODZ(IMAG).GT.0.0D0) THEROT(IMAG)=-THEROT(IMAG)
          ENDDO   !IMAG=1,NMAGMOD
        ENDIF !ITHESYMD

        DO IMAG=1,NMAGMOD

          RADIMOD(IMAG)=RADIMOD(IMAG)*SCALRAD
          THEROT(IMAG)=THEROT(IMAG)*SCALTHE

          IF (CENMODY(IMAG).LT.0.0D0) THEN
            IF (CENMODX(IMAG).GE.0.0D0) THEN
              THEROT(IMAG)=THEROT(IMAG)+ITHESYML*ITHESYMD*THEGROTL
            ELSE
              THEROT(IMAG)=THEROT(IMAG)+ITHESYML*THEGROTL
            ENDIF
          ELSE
            IF (CENMODX(IMAG).GE.0.0D0) THEN
              THEROT(IMAG)=THEROT(IMAG)+ITHESYMD*THEGROTU
            ELSE
              THEROT(IMAG)=THEROT(IMAG)+THEGROTL
            ENDIF
          ENDIF

        ENDDO  !IMAG=1,NMAGMOD

        ICAL=1

      ENDIF !(ICAL.EQ.0)

      DO IMAG=1,NMAGMOD

        BCCOS=BCMOD(IMAG)*DCOSD(THEROT(IMAG))*SCALMOD
        IF (ABS(BCCOS).LT.1.0D-10) BCCOS=0.0D0
        BCSIN=BCMOD(IMAG)*DSIND(THEROT(IMAG))*SCALMOD
        IF (ABS(BCSIN).LT.1.0D-10) BCSIN=0.0D0

        RHO2=RADIMOD(IMAG)*RADIMOD(IMAG)
        DLY=2.D0*RADIMOD(IMAG)/NSLICE

        DO I=1,NSLICE

          KMAG=IMAG+I-1
          LMAG=KMAG+NSLICE

          Y2=CENMODY(IMAG)-RADIMOD(IMAG)+DLY*I
          Y1=Y2-DLY
          Y=(Y2+Y1)/2.D0

          XMOD(KMAG)=CENMODX(IMAG)
          DXMOD(KMAG)=2.D0*DSQRT(RHO2-(Y-CENMODY(IMAG))**2)
          YMOD(KMAG)=Y
          DYMOD(KMAG)=DLY
          ZMOD(KMAG)=CENMODZ(IMAG)
          DZMOD(KMAG)=ZLENMOD(IMAG)
          BCBCMOD(KMAG)=BCCOS
          THEMOD(KMAG)=0.D0
          PHIMOD(KMAG)=0.D0

          XMOD(LMAG)=XMOD(KMAG)
          YMOD(LMAG)=YMOD(KMAG)
          ZMOD(LMAG)=ZMOD(KMAG)
          DXMOD(LMAG)=DXMOD(KMAG)
          DYMOD(LMAG)=DYMOD(KMAG)
          DZMOD(LMAG)=DZMOD(KMAG)
          BCBCMOD(LMAG)=BCSIN
          THEMOD(LMAG)=90.D0
          PHIMOD(LMAG)=0.D0

          IMAGTOT=IMAGTOT+1
          DX(IMAGTOT)=XMOD(KMAG)
          DY(IMAGTOT)=YMOD(KMAG)
          DZ(IMAGTOT)=ZMOD(KMAG)
          XLEN(IMAGTOT)=DXMOD(KMAG)
          YLEN(IMAGTOT)=DYMOD(KMAG)
          ZLEN(IMAGTOT)=DZMOD(KMAG)
          THETA(IMAGTOT)=THEMOD(KMAG)*GRARAD
          PHI(IMAGTOT)=PHIMOD(KMAG)*GRARAD
          BC(IMAGTOT)=BCBCMOD(KMAG)

          IMAGTOT=IMAGTOT+1
          DX(IMAGTOT)=XMOD(LMAG)
          DY(IMAGTOT)=YMOD(LMAG)
          DZ(IMAGTOT)=ZMOD(LMAG)
          XLEN(IMAGTOT)=DXMOD(LMAG)
          YLEN(IMAGTOT)=DYMOD(LMAG)
          ZLEN(IMAGTOT)=DZMOD(LMAG)
          THETA(IMAGTOT)=THEMOD(LMAG)*GRARAD
          PHI(IMAGTOT)=PHIMOD(LMAG)*GRARAD
          BC(IMAGTOT)=BCBCMOD(LMAG)

        ENDDO   !NSLICE

      ENDDO !NMAGMOD

      RETURN
      END
