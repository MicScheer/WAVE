*CMZ :  2.67/02 08/05/2012  12.08.50  by  Michael Scheer
*CMZ :  1.01/00 07/10/2004  09.04.51  by  Michael Scheer
*CMZ :  1.00/03 04/10/2004  15.59.52  by  Michael Scheer
*-- Author :    Michael Scheer   01/10/2004
      SUBROUTINE MSHSYNCSPEC(MODE,Y,PSI,WRAD_C)
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

C--- CALCULATE normalized DIPOL SPECTRUM (G1)

      IMPLICIT NONE

      EXTERNAL FUNCTION DBSKR3
      DOUBLE PRECISION Y,PSI,WRAD_C,PER,PAR,XI,XX,XX1,BK13,BK23,DBSKR3

      INTEGER MODE

      IF (Y.LE.0.D0) THEN
        STOP '*** Error:     E/Ec is zero in MSHSYNCSPEC'
      ENDIF

      IF (MODE.EQ.1) THEN

        XX=PSI**2
        XX1=XX+1.D0
        XI=Y*DSQRT(XX1)**3/2.D0
        BK13=DBSKR3(XI,1)
        BK23=DBSKR3(XI,2)

        PAR=(Y*XX1*BK23)**2
        PER=Y*Y*XX*XX1*BK13**2

        WRAD_C=PAR+PER

      ELSE IF (MODE.EQ.-1) THEN

C USE ANALYTICAL FIT TO G1

        IF(Y.LT.4.) THEN
          WRAD_C = 391.8 * Y**0.333 * EXP(-Y*0.8307)-
     &      192.0 * Y**0.500 * EXP(-Y*0.7880)
        ELSE
          WRAD_C = 164.0 * Y**0.500 * EXP(-Y)
        END IF

        wrad_c=wrad_c * 0.007565 ! anders als in geant3 wird hier G1 berechnet

      ELSE !MODE

        STOP '*** ERROR IN MSHSYNCSPEC'

      ENDIF !MODE

      RETURN
      END
