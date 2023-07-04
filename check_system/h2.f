*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.33/09 09/05/2001  16.51.33  by  Michael Scheer
*CMZ :  2.33/00 02/05/2001  14.43.10  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.00/00 15/12/98  10.35.46  by  Michael Scheer
*CMZ : 00.02/02 15/01/97  15.18.41  by  Michael Scheer
*CMZ : 00.00/07 18/05/94  14.54.15  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.54.12  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.24  by  Michael Scheer
*-- Author : Michael Scheer
      DOUBLE PRECISION FUNCTION H2(Y)
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

C--- FUNCTION H2(Y)

      IMPLICIT NONE

      DOUBLE PRECISION Y,XI,BK23,DBSKR3

      IF (Y.EQ.0.D0) THEN
          WRITE(6,*)'     E/Ec is zero in function H2'
          STOP '--- PROGRAM WAVE ABORTED ---'
      ENDIF

       XI=Y/2.D0
       BK23=DBSKR3(XI,2)
       H2=(Y*BK23)**2

      RETURN
      END
