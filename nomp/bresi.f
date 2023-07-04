*CMZ : 00.01/09 05/10/95  16.28.55  by  Michael Scheer
*CMZ :  1.00/03 27/09/95  16.41.21  by  Michael Scheer
*CMZ :  1.00/02 26/09/95  17.45.08  by  Michael Scheer
*CMZ :  1.00/01 22/09/95  18.24.37  by  Michael Scheer
*-- Author :    Michael Scheer   22/09/95

      SUBROUTINE BRESI(NPOI,NPOIX,NPOIY,NPOIZ
     &                  ,X,Y,Z,BX,BY,BZ,BXF,BYF,BZF
     &                  ,RESBX,RESBY,RESBZ,RESB
     &                  ,BXAMEAN,BYAMEAN,BZAMEAN,BAMEAN
     &                  ,BERRMX,BXERRMX,BYERRMX,BZERRMX)
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

      IMPLICIT NONE

      INTEGER NPOI,I,NPOIX,NPOIY,NPOIZ

      DOUBLE PRECISION
     &                  BX(NPOI),BY(NPOI),BZ(NPOI)
     &                  ,X(NPOI),Y(NPOI),Z(NPOI)
     &                  ,BXF(NPOI),BYF(NPOI),BZF(NPOI)
     &                  ,RESBX,RESBY,RESBZ,RESB
     &                  ,BXAMEAN,BYAMEAN,BZAMEAN,BAMEAN
     &                  ,BBX,BBY,BBZ,BBXF,BBYF,BBZF,BB,BBF
     &                  ,DBX,DBY,DBZ,DBX2,DBY2,DBZ2,DBB
     &                  ,BERRMX(8),BXERRMX(8),BYERRMX(8),BZERRMX(8)

      RESBX=0.D0
      RESBY=0.D0
      RESBZ=0.D0
      RESB=0.D0

      BXAMEAN=0.D0
      BYAMEAN=0.D0
      BZAMEAN=0.D0
      BAMEAN=0.D0

      BERRMX(8)=0.D0
      BXERRMX(8)=0.D0
      BYERRMX(8)=0.D0
      BZERRMX(8)=0.D0

      DO I=1,NPOI

          BBX=BX(I)
          BBY=BY(I)
          BBZ=BZ(I)
          BBXF=BXF(I)
          BBYF=BYF(I)
          BBZF=BZF(I)

          IF (BBX.EQ.-9999.) THEN
         BBX=0.D0
         BBXF=0.D0
          ENDIF

          IF (BBY.EQ.-9999.) THEN
         BBY=0.D0
         BBYF=0.D0
          ENDIF

          IF (BBZ.EQ.-9999.) THEN
         BBZ=0.D0
         BBZF=0.D0
          ENDIF

          BB=DSQRT(BBX*BBX+BBY*BBY+BBZ*BBZ)
          BBF=DSQRT(BBXF*BBXF+BBYF*BBYF+BBZF*BBZF)

          DBX=DABS(BBX-BBXF)
          DBY=DABS(BBY-BBYF)
          DBZ=DABS(BBZ-BBZF)
          DBB=DABS(BB-BBF)

          DBX2=DBX**2
          DBY2=DBY**2
          DBZ2=DBZ**2

          IF (DBB.GT.DABS(BERRMX(8))) THEN
         BERRMX(1)=X(I)
         BERRMX(2)=Y(I)
         BERRMX(3)=Z(I)
         BERRMX(4)=BX(I)
         BERRMX(5)=BY(I)
         BERRMX(6)=BZ(I)
         BERRMX(7)=BB
         BERRMX(8)=BB-BBF
          ENDIF
          IF (DBX.GT.DABS(BXERRMX(8))) THEN
         BXERRMX(1)=X(I)
         BXERRMX(2)=Y(I)
         BXERRMX(3)=Z(I)
         BXERRMX(4)=BX(I)
         BXERRMX(5)=BY(I)
         BXERRMX(6)=BZ(I)
         BXERRMX(7)=BB
         BXERRMX(8)=BBX-BBXF
          ENDIF
          IF (DBY.GT.DABS(BYERRMX(8))) THEN
         BYERRMX(1)=X(I)
         BYERRMX(2)=Y(I)
         BYERRMX(3)=Z(I)
         BYERRMX(4)=BX(I)
         BYERRMX(5)=BY(I)
         BYERRMX(6)=BZ(I)
         BYERRMX(7)=BB
         BYERRMX(8)=BBY-BBYF
          ENDIF
          IF (DBZ.GT.DABS(BZERRMX(8))) THEN
         BZERRMX(1)=X(I)
         BZERRMX(2)=Y(I)
         BZERRMX(3)=Z(I)
         BZERRMX(4)=BX(I)
         BZERRMX(5)=BY(I)
         BZERRMX(6)=BZ(I)
         BZERRMX(7)=BB
         BZERRMX(8)=BBZ-BBZF
          ENDIF

          BXAMEAN=BXAMEAN+DABS(BBX)
          BYAMEAN=BYAMEAN+DABS(BBY)
          BZAMEAN=BZAMEAN+DABS(BBZ)
          BAMEAN=BAMEAN+DSQRT(BBX*BBX+BBY*BBY+BBZ*BBZ)

          RESBX=RESBX+DBX2
          RESBY=RESBY+DBY2
          RESBZ=RESBZ+DBZ2
          RESB=RESB
     &     +((DSQRT(BBX**2+BBY**2+BBZ**2)
     &      -DSQRT(BBXF**2+BBYF**2+BBZF**2)))**2

      ENDDO !NPOI

      BAMEAN=BAMEAN/NPOI
      BXAMEAN=BXAMEAN/NPOIX
      BYAMEAN=BYAMEAN/NPOIY
      BZAMEAN=BZAMEAN/NPOIZ

      RESB=DSQRT(RESB/NPOI)
      IF (NPOIX.NE.0) RESBX=DSQRT(RESBX/NPOIX)
      IF (NPOIY.NE.0) RESBY=DSQRT(RESBY/NPOIY)
      IF (NPOIZ.NE.0) RESBZ=DSQRT(RESBZ/NPOIZ)

      RETURN
      END
