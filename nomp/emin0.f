*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  16.22.05  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.50.01  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.44  by  Michael Scheer
*-- Author : Michael Scheer
C****************************************************************
      SUBROUTINE EMIN0(CD,C5,C2,DI2,DI5,ENERG,EMINEM,EMIN,IOK)
C****************************************************************
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
C     LOEST GLEICHUNGSSYSTEM FUER MINIMALES E
          IMPLICIT NONE
          DOUBLE PRECISION CD,C5,C2,DI2,DI5,E,DEFUN,EFUN,EMINEM,EMIN,ENERG
          DOUBLE PRECISION A,B,EN
          INTEGER I,IOK

      ENERG=ENERG

C--- NEWTON-VERFAHREN

      IOK=1
      E=EMINEM
      DO I=1,100
          A=DEFUN(CD,C5,C2,DI2,DI5,E)
          B=EFUN(CD,C5,C2,DI2,DI5,E)-A*E
          EN=-B/A
          IF(DABS((EN-E)/EN).LT.1.D-10) GOTO 100
          E=EN
          IF(E.GT.1.D6) THEN  !FAENGT OVERFLOW AB
         IOK=0
              RETURN
          ENDIF
          END DO
      IOK=0
      RETURN
C     STOP '*** S/R EMIN0: KEINE LOESUNG GEFUNDEN ***'
100   CONTINUE
      EMIN=EN
      IF(EMIN.LT.0.) IOK=0
      RETURN
      END
