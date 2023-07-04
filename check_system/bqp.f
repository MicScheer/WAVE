*CMZ :  2.67/01 15/03/2012  10.51.33  by  Michael Scheer
*CMZ :  2.66/07 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  1.02/00 18/12/97  13.35.57  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.24.44  by  Michael Scheer
*CMZ : 00.00/05 29/04/94  19.35.51  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.43  by  Michael Scheer
*-- Author : Michael Scheer
C***********************************************************************
      SUBROUTINE BQP(XI,YI,ZI,BX,BY,BZ,IMAG)
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
      INTEGER IMAG

*KEEP,mgsqc.
      include 'mgsqc.cmn'
*KEND.

      DOUBLE PRECISION Y,Z,BX,BY,BZ,G,XI,YI,ZI,PHI,SPHI,CPHI,XCEN,ZCEN,DX,DZ
      DOUBLE PRECISION DXA,DZA,DXE,DZE,XA,XE,ZA,ZE,ADUMA,ADUME

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      PHI=PMAG(5,IMAG)
      CPHI=DCOS(PHI)
      SPHI=DSIN(PHI)
      XCEN=PMAG(3,IMAG)
      ZCEN=PMAG(4,IMAG)

      XA=XCEN-CPHI*PMAG(1,IMAG)/2.
      ZA=ZCEN-SPHI*PMAG(1,IMAG)/2.
      XE=XCEN+CPHI*PMAG(1,IMAG)/2.
      ZE=ZCEN+SPHI*PMAG(1,IMAG)/2.
      DXA=XI-XA
      DZA=ZI-ZA
      DXE=XI-XE
      DZE=ZI-ZE

      ADUMA=DXA*CPHI+DZA*SPHI
      ADUME=DXE*CPHI+DZE*SPHI

      IF (IWFILF.NE.99.AND.ADUMA*ADUME.LE.0.) THEN

        DX=XI-XCEN
        DZ=ZI-ZCEN
        Z=-SPHI*DX+CPHI*DZ
        Y=YI

        G=PMAG(2,IMAG)*EMOM/CLIGHT1

        BX=0.
        BY=G*Z
        BZ=G*Y

      ELSE

        BX=0.
        BY=0.
        BZ=0.

      ENDIF !IWFILF

      RETURN
      END
