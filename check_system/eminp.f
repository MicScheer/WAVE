*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.15/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.36.44  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.45  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE EMINP(R,B3,T0,TK,EMIN)
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

C     SOLVES EQUATION TO DETERMINE MINIMUM ENERGY FOR GIVEN POLARIZATION LEVEL
C     LOGBUCH SEITE 189

          IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

          DOUBLE PRECISION R,B3,T0,TK,EMIN,E
          DOUBLE PRECISION A,B,EN,DEFUNP,EFUNP
          INTEGER I

C--- ALGORITHM OF NEWTON X(N+1)=X(N)-F(X(N))/F'(X(N))

      E=1.D0
      DO I=1,100
          A=DEFUNP(R,B3,T0,TK,E)
          B=EFUNP(R,B3,T0,TK,E)-A*E
          EN=-B/A
          IF(DABS((EN-E)/EN).LT.1.D-10) GOTO 100
          E=EN
          END DO

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'*** WARNING SR EMINP ***'
      WRITE(LUNGFO,*)'NEWTON ALGORITHM FAILED TO DETERMINE MINIMUM ENERGY'
      WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

100   CONTINUE

      EMIN=EN

      RETURN
      END
