*CMZ :  2.41/10 14/08/2002  17.34.02  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  1.03/06 06/08/98  18.01.12  by  Michael Scheer
*CMZ :  1.00/00 03/06/97  13.54.04  by  Michael Scheer
*CMZ : 00.00/04 07/09/95  10.59.03  by  Michael Scheer
*-- Author :    Michael Scheer   05/09/95

      SUBROUTINE SPLITMAGCYL
     &  (ZLEN,BC,RHO,CENX,CENY,DZ0,THEROT,NCYL,LUNMAG,FILENAME)
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

      EXTERNAL DCOSD,DSIND
      DOUBLE PRECISION DCOSD,DSIND

      CHARACTER(72) FILENAME

      INTEGER I,NCYL,LUNMAG

      DOUBLE PRECISION RHO,CENY,CENX,THETA,Y0,Y,RHO2,BCSIN,BCCOS
      DOUBLE PRECISION THEROT,PHI,BC,XLEN,YLEN,ZLEN,DX0,DY0,DZ0
      DOUBLE PRECISION DX,X1,X2,DY,Y1,Y2

      DATA THETA/0.0D0/
      DATA PHI/0.0D0/

      BCCOS=BC*DCOSD(THEROT)
      BCSIN=BC*DSIND(THEROT)

      OPEN(UNIT=LUNMAG,FILE=FILENAME,STATUS='NEW')

      WRITE(LUNMAG,*)2*NCYL

      RHO2=RHO*RHO
      DY=2.D0*RHO/NCYL
      Y0=CENY-RHO
      DX0=CENX
      YLEN=DY
      DO I=1,NCYL
          Y2=Y0+DY*I
          Y1=Y2-DY
          Y=(Y2+Y1)/2.D0
          DX=DSQRT(RHO2-(Y-CENY)**2)
          X1=CENX-DX
          X2=CENX+DX
          XLEN=2.D0*DX
          DY0=Y
          WRITE(LUNMAG,1000)XLEN,YLEN,ZLEN
          WRITE(LUNMAG,1000)THETA,PHI,BCCOS
          WRITE(LUNMAG,1000)DX0,DY0,DZ0
          WRITE(LUNMAG,1000)XLEN,YLEN,ZLEN
          WRITE(LUNMAG,1000)THETA+90.D0,PHI,BCSIN
          WRITE(LUNMAG,1000)DX0,DY0,DZ0
      ENDDO !NCYL

      CLOSE(LUNMAG)

      RETURN

1000  FORMAT (3(1PD25.13))

      END
