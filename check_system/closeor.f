*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.49.07  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.29  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE CLOSEOR(A01,A10,A11,A20,A02,X,P)
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

      DOUBLE PRECISION A01,A10,A11,A20,A02,X,P

      X=0.
      P=0.

      IF(A11.NE.1.D0) THEN
          X=(A01*(1.-A11)+2.*A02*A10)/((1.-A11)**2-4.*A02*A20)
          P=A10/(1.-A11)+2.*A20/(1.-A11)*X
C     ELSEIF(A11.EQ.1.D0)THEN
C         IF (A20.NE.0) X=-A10/(2.*A20)
C         IF (A02.NE.0) P=-A01/(2.*A02)
C         IF (A02.EQ.0. .AND. A01.NE.0. .OR. A20.EQ.0. .AND. A10.NE.0.)
C     &      STOP '*** S/R CLOSEOR: KOEFFIZIENTEN INKONSISTENT ***'
      ENDIF

c     P=P/DSQRT(1.D0-P*P)

      RETURN
      END
