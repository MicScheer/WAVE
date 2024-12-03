*CMZ :          18/11/2024  10.51.09  by  Michael Scheer
*CMZ :  3.04/00 05/01/2018  16.13.44  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.65/03 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.63/02 25/03/2008  09.35.21  by  Michael Scheer
*CMZ :  2.52/11 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  1.02/00 19/12/97  17.58.18  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  14.07.41  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.46.47  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.43  by  Michael Scheer
*-- Author : Michael Scheer
C***********************************************************************
      SUBROUTINE BDI(XI,YI,ZI,BX,BY,BZ,IMAGI)

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

      INTEGER IMAG,imagi

*KEEP,mgsqc.
      include 'mgsqc.cmn'
*KEND.

      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,BY0,BY1,BY2,XLEN2,XI,YI,ZI,AY1,AY2

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

c pmag(1,imag): deflection angle (rad)
c pmag(2,imag): bending radius (T)
c pmag(3,imag): Center of magnet (m)
c pmag(4,imag): Width of edge

c 'DI' type (imagi>0)
c pmag(5,imag): cos(rotation angle)
c pmag(6,imag): sin(rotation angle)

c 'DIF' type (imagi<0)
c pmag(7,imag): cos(rotation angle)
c pmag(8,imag): sin(rotation angle)

      imag=iabs(imagi)

      IF (
     &    (IWFILF.EQ.99. .OR. IMGSQF.EQ.0)
     &    .AND.
     &    PMAG(1,IMAG)*PMAG(2,IMAG).NE.0.
     &    ) THEN

        if (imagi.gt.0) then
          X=pmag(5,imag)*(XI-PMAG(3,IMAG))+pmag(6,imag)*zi
          Y=YI
          z=-pmag(6,imag)*(XI-PMAG(3,IMAG))+pmag(5,imag)*zi
        else
          X=pmag(7,imag)*(XI-PMAG(3,IMAG))+pmag(7,imag)*zi
          Y=YI
          z=-pmag(8,imag)*(XI-PMAG(3,IMAG))+pmag(8,imag)*zi
        endif

        BX=0.
        BZ=0.

        BY0=EMOM/CLIGHT1/PMAG(2,IMAG)
        XLEN2=DABS(PMAG(2,IMAG)*sin(pmag(1,imag)/2.0d0))


        AY1=(+X-XLEN2)*PMAG(4,IMAG)
        AY2=(-X-XLEN2)*PMAG(4,IMAG)

        IF (AY1.GT.70.0D0) THEN
          BY1=0.0d0
        ELSE IF (AY1.LT.-70.) THEN
          BY1=1.0D0
        ELSE
          BY1=1.0D0/(1.0D0+DEXP(AY1))
        ENDIF

        IF (AY2.GT.70.0D0) THEN
          BY2=0.0D0
        ELSE IF (AY2.LT.-70.0D0) THEN
          BY2=1.0D0
        ELSE
          BY2=1.0D0/(1.0D0+DEXP(AY2))
        ENDIF

        BY=BY0*BY1*BY2*CORR(IMAG)

      ELSE

        BX=0.0d0
        BY=0.0d0
        BZ=0.0d0

      ENDIF

      RETURN
      END
