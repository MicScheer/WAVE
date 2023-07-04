*CMZ :  3.03/00 13/07/2015  15.41.18  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.34/09 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.17/00 03/11/2000  10.22.05  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  18.04.39  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.36  by  Michael Scheer
*CMZ :  2.13/03 16/12/99  10.48.38  by  Michael Scheer
*CMZ :  1.03/06 11/06/98  18.24.32  by  Michael Scheer
*CMZ : 00.01/04 28/11/94  18.44.59  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.57.19  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.20  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WPOLY2(ISOUR,IFREQ)
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
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

C--- DETERMINE 2D CUBIC POLYNOMIALS WHICH DESCRIBE THE INTENSITY INSIDE A
C    MASH OF THE PINHOLE
C    THE USED SR WBCUCOF STEMS FORM NUMERICAL RECIPIES (BCUCOF)
C    THE POLYNOMIAL COEFFICIENTS OF THE MASHES ARE STORED

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEND.

      INTEGER ISOUR,IFREQ,IOBSV,IZ,IY,NDIA,NPL,IPL,IDIA,I1,I2,I3,I4
      INTEGER IYY,IZZ,ISS,NOBSVZ1

      DOUBLE PRECISION YPP1(4),YPP2(4),YPP3(4),CSP,DSP,YPP0,YPPN

      REAL*4 Y(4),Y1(4),Y2(4),Y12(4),C(4,4)

      ALLOCATE(IOBUFF(NOBSV))

C- SPLINES IN Z

      DO IZ=1,MAX(NOBSVZ,NOBSVY)
          DOBUFF1(IZ)=DFLOAT(IZ-1)
      ENDDO !IZ

        IOBSV=0
        DO IY=1,NOBSVY

        DO IZ=1,NOBSVZ

          IOBSV=IOBSV+1
          DOBUFF(IZ)=SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
        ENDDO      !IZ

           YPP0=1.D30
           YPPN=1.D30
C       CALL FSPLINEZ(DOBUFF1,DOBUFF,NOBSVZ,YPP0,YPPN,DOBUFF2)
C060793    CALL FSPLINDX(DOBUFF1,            S,NOBSVZ,YPP0,YPPN,DOBUFF2)
        CALL FSPLINDX(DOBUFF1(2)-DOBUFF1(1),DOBUFF,NOBSVZ,0.D0,0.D0,DOBUFF2)

        DO IZ=1,NOBSVZ

          SPCOEFU(1,IZ+(IY-1)*NOBSVZ)=DOBUFF2(IZ)


        ENDDO      !IZ
        ENDDO      !IY

C- SPLINES IN Y


        DO IZ=1,NOBSVZ
        DO IY=1,NOBSVY

          IOBSV=(IY-1)*NOBSVZ+IZ
          DOBUFF(IY)=SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))

        ENDDO      !IY

           YPP0=1.D30
           YPPN=1.D30
C       CALL FSPLINEY(DOBUFF1,DOBUFF,NOBSVY,YPP0,YPPN,DOBUFF2)
C060793    CALL FSPLINDX(DOBUFF1,            S,NOBSVY,YPP0,YPPN,DOBUFF2)
        CALL FSPLINDX(DOBUFF1(2)-DOBUFF1(1),DOBUFF,NOBSVY,0.D0,0.D0,DOBUFF2)

        DO IY=1,NOBSVY

          SPCOEFU(2,IZ+(IY-1)*NOBSVZ)=DOBUFF2(IY)

        ENDDO      !IY
        ENDDO      !IZ

C- SPLINES IN YZ-DIRECTION

      DO IY=1,MAX(NOBSVY,NOBSVZ)
        DOBUFF1(IY)=DSQRT(2.D0)*DFLOAT(IY-1)
      ENDDO !IY

        NDIA=NOBSVZ+NOBSVY-1  !NUMBER OF LINES ACROSS THE PINHOLE

        DO IDIA=1,NDIA

C- NUMBER OF POINTS OF CURRENT LINE

          IF (IDIA.LE.NOBSVY) THEN
         NPL=MIN(IDIA,NOBSVZ)
         DO IPL=1,NPL
             IOBSV=NOBSVY*NOBSVZ-NOBSVZ*IDIA+(IPL-1)*NOBSVZ+IPL
             IOBUFF(IPL)=IOBSV
              DOBUFF(IPL)=SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
              ENDDO   !IPL
          ELSE
         NPL=NDIA-IDIA+1

CORR 2.10.92
         NPL=MIN(NPL,NOBSVY)
CORR 2.10.92

         DO IPL=1,NPL
             IOBSV=IDIA-NOBSVY+(IPL-1)*NOBSVZ+IPL
             IOBUFF(IPL)=IOBSV
              DOBUFF(IPL)=SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
              ENDDO   !IPL
          ENDIF

           YPP0=1.D30
           YPPN=1.D30
          IF (NPL.LE.NDOBSVZ) THEN
C            CALL FSPLINEZ(DOBUFF1,DOBUFF,NPL,YPP0,YPPN,DOBUFF2)
C060793         CALL FSPLINDX(DOBUFF1,            S,NPL,YPP0,YPPN,DOBUFF2)
                  CALL FSPLINDX(DOBUFF1(2)-DOBUFF1(1),DOBUFF,NPL,0.D0,0.D0,DOBUFF2)
          ELSEIF (NPL.LE.NDOBSVY) THEN
C            CALL FSPLINEY(DOBUFF1,DOBUFF,NPL,YPP0,YPPN,DOBUFF2)
C060793         CALL FSPLINDX(DOBUFF1,            S,NPL,YPP0,YPPN,DOBUFF2)
                  CALL FSPLINDX(DOBUFF1(2)-DOBUFF1(1),DOBUFF,NPL,0.D0,0.D0,DOBUFF2)
          ELSE
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** ERROR IN WPOLY2 ***'
             WRITE(LUNGFO,*)
     & 'SEVERE ERROR! CHECK SOURCE CODE USING THE DEBUGGER, SORRY'
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)'*** ERROR IN WPOLY2 ***'
             WRITE(6,*)
     & 'SEVERE ERROR! CHECK SOURCE CODE USING THE DEBUGGER, SORRY'
             WRITE(6,*)
             STOP
          ENDIF

          DO IPL=1,NPL

            ISS=IOBUFF(IPL)
            SPCOEFU(3,ISS)=DOBUFF2(IPL)

          ENDDO   !IPL

          ENDDO   !IDIA

C--- NOW THE COEFFICIENTS OF THE BICUBIC SPLINES (NUMERICAL RECIPIES PAGE 98)


C- LOOP OVER MASHES

        IF (IF1DIM.EQ.0) THEN
          NOBSVZ1=NOBSVZ-1
        ELSE
          NOBSVZ1=1
        ENDIF   !IF1DIM

          DO IY=1,NOBSVY-1
          DO IZ=1,NOBSVZ1

             IF (IF1DIM.EQ.0) THEN
                I1=IZ+NOBSVZ*(IY-1)
                I2=I1+1
                I4=I1+NOBSVZ
                I3=I4+1
             ELSE
                I1=IZ+NOBSVZ*(IY-1)
                I2=I1
                I4=I1+NOBSVZ
                I3=I4
             ENDIF   !IF1DIM

             Y(1)=SPEC(ISOUR+NSOURCE*(I1-1+NOBSV*(IFREQ-1)))    !INTENSITIES AT MASH POINTS
             Y(2)=SPEC(ISOUR+NSOURCE*(I2-1+NOBSV*(IFREQ-1)))    !INTENSITIES AT MASH POINTS
             Y(3)=SPEC(ISOUR+NSOURCE*(I3-1+NOBSV*(IFREQ-1)))    !INTENSITIES AT MASH POINTS
             Y(4)=SPEC(ISOUR+NSOURCE*(I4-1+NOBSV*(IFREQ-1)))    !INTENSITIES AT MASH POINTS
             YPP1(1)=SPCOEFU(1,I1)    !SECOND DERIVATIVES WITH
             YPP1(2)=SPCOEFU(1,I2)    !RESPECT TO Z
             YPP1(3)=SPCOEFU(1,I3)
             YPP1(4)=SPCOEFU(1,I4)

             CSP=-2.D0/6.D0
             DSP=-1.D0/6.D0
             Y1(1)=-Y(1)+Y(2)+CSP*YPP1(1)+DSP*YPP1(2)
             Y1(4)=-Y(4)+Y(3)+CSP*YPP1(4)+DSP*YPP1(3)

             CSP=+1.D0/6.D0
             DSP=+2.D0/6.D0
             Y1(2)=-Y(1)+Y(2)+CSP*YPP1(1)+DSP*YPP1(2)
             Y1(3)=-Y(4)+Y(3)+CSP*YPP1(4)+DSP*YPP1(3)

             YPP2(1)=SPCOEFU(2,I1)    !SECOND DERIVATIVES WITH
             YPP2(2)=SPCOEFU(2,I2)    !RESPECT TO Y
             YPP2(3)=SPCOEFU(2,I3)
             YPP2(4)=SPCOEFU(2,I4)


             CSP=-2.D0/6.D0
             DSP=-1.D0/6.D0
             Y2(1)=-Y(1)+Y(4)+CSP*YPP2(1)+DSP*YPP2(4)
             Y2(2)=-Y(2)+Y(3)+CSP*YPP2(2)+DSP*YPP2(3)

             CSP=+1.D0/6.D0
             DSP=+2.D0/6.D0
             Y2(4)=-Y(1)+Y(4)+CSP*YPP2(1)+DSP*YPP2(4)
             Y2(3)=-Y(2)+Y(3)+CSP*YPP2(2)+DSP*YPP2(3)

             YPP3(1)=SPCOEFU(3,I1)    !DERIVATIVES WITH RESPECT
             YPP3(2)=SPCOEFU(3,I2)    !TO Z AND Y
             YPP3(3)=SPCOEFU(3,I3)
             YPP3(4)=SPCOEFU(3,I4)

             IF(IF1DIM.EQ.0) THEN
                Y12(1)=YPP3(1)-(YPP1(1)+YPP2(1))/2.D0
                Y12(2)=YPP3(2)-(YPP1(2)+YPP2(2))/2.D0
                Y12(3)=YPP3(3)-(YPP1(3)+YPP2(3))/2.D0
                Y12(4)=YPP3(4)-(YPP1(4)+YPP2(4))/2.D0
             ELSE
                Y12(1)=0.
                Y12(2)=0.
                Y12(3)=0.
                Y12(4)=0.
             ENDIF   !IF1DIM

             CALL WBCUCOF(Y,Y1,Y2,Y12,1.E0,1.E0,C)

             DO IYY=1,4
             DO IZZ=1,4
           COFOLD(IZZ,IYY,I1)=C(IZZ,IYY)
             ENDDO   !IZZ
             ENDDO   !IYY


          ENDDO   !IY
          ENDDO   !IZ

      DEALLOCATE(IOBUFF)
      RETURN
      END
