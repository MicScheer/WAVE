*CMZ :  4.00/11 19/04/2021  16.49.53  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.54/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.41.44  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.01/08 31/05/95  13.30.40  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  16.00.44  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  18.02.31  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.57  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BTABxyz(XIN,DUMY,DUMZ,BX,BY,BZ,AX,AY,AZ)
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

C     BTABxyz reads By and Bz from data files and interpolatEs fields

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      INTEGER ICAL

      DOUBLE PRECISION XIN,DUMY,DUMZ,DUMZY,XS,XE,XS1,XS2,XE1,XE2
      DOUBLE PRECISION AX,AY,AZ,BX,BY,BZ
      DOUBLE PRECISION AXX,AYY,AZZ,BXX,BYY,BZZ
      DOUBLE PRECISION AXXX,AYYY,AZZZ,BXXX,BYYY,BZZZ

      DATA ICAL/0/
      DATA XS1,XS2,XE1,XE2/4*0.0/

      DUMZY=DUMY
      DUMZY=DUMZ

      BXX=0.
      BYY=0.
      BZZ=0.

      AXX=0.
      AYY=0.
      AZZ=0.

      BXXX=0.
      BYYY=0.
      BZZZ=0.

      AXXX=0.
      AYYY=0.
      AZZZ=0.

      XS=XSTART
      XE=XSTOP

      CALL BTAB (XIN,0.D0,0.D0, BXX, BYY, BZZ, AXX, AYY, AZZ)

      IF (ICAL.EQ.0) THEN

         XS1=XSTART
         XE1=XSTOP

      ENDIF !ICAL

      CALL BTABX(XIN,0.D0,0.D0,BX,BY,BZ,AX,AY,AZ,XS,XE)

      IF (ICAL.EQ.0) THEN

          XS2=XSTART
          XE2=XSTOP

          IF (XS1.NE.XS2) XSTART=DMIN1(XS1,XS2)
          IF (XE1.NE.XE2)  XSTOP=DMAX1(XE1,XE2)

      ENDIF !ICAL

      CALL BTABZ(XIN,0.D0,0.D0,BXXX,BYYY,BZZZ,AXXX,AYYY,AZZZ,XS,XE)

      IF (ICAL.EQ.0) THEN

          XS2=XSTART
          XE2=XSTOP

          IF (XS1.NE.XS2) XSTART=DMIN1(XS1,XS2)
          IF (XE1.NE.XE2)  XSTOP=DMAX1(XE1,XE2)

          ICAL=1

      ENDIF !ICAL

      BX=bx+BXX+BXXX
      BY=by+BYY+BYYY
      BZ=bz+BZZ+BZZZ

      AX=ax+AXX+AXXX
      AY=ay+AYY+AYYY
      AZ=az+AZZ+AZZZ

      RETURN
      END
