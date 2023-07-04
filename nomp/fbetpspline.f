*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.41.44  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  16.36.10  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.50.53  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.22  by  Michael Scheer
*-- Author : Michael Scheer
C******************************************************************
      SUBROUTINE FBETPSPLINE(XI,YO)
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

C--- DOUBLE PRECISION SUBROUTINE FOR CUBIC SPLINE INTERPOLATION
C    SUBROUTINES FROM NUMERICAL RECIEPIES

      IMPLICIT NONE
      INTEGER NPOINT,ICAL,NMAX
      PARAMETER (NMAX=100) !BEACHTE AUCH NMAX IN S/R SPLINE UND SPLINT
      DOUBLE PRECISION X(NMAX),Y(NMAX),Y2(NMAX),XI,YO

      DATA ICAL/0/

C--- DATA SET I.E. TABULATED FUNCTION Y(X)

      DATA NPOINT   /       6        /

      DATA X(1),Y(1)/       1.0,    0.4704       /
      DATA X(2),Y(2)/       2.0,    0.3332       /
      DATA X(3),Y(3)/       3.0,    0.2405       /
      DATA X(4),Y(4)/       4.0,    0.2018       /
      DATA X(5),Y(5)/       6.0,    0.1329       /
      DATA X(6),Y(6)/       8.0,    0.08906      /

      IF (ICAL.NE.1) THEN
          CALL SPLINE(X,Y,NPOINT,2.D30,2.D30,Y2)
          ICAL=1
      ENDIF

      CALL SPLINT(X,Y,Y2,NPOINT,XI,YO)

      RETURN
      END
