*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.63/02 05/02/2008  14.00.22  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.52  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE PARABEL(x1,y1,x2,y2,x3,y3,A)

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
C--- CALCULATES A(1), A(2), A(3) OF PARABOLA Y=A1+A2*X+A3*X**2 FROM COORDINATES
C    OF THREE POINTS

      IMPLICIT NONE

        DOUBLE PRECISION x1,x2,x3,y1,y2,y3,a(3)
        DOUBLE PRECISION x(3),y(3),yp(3),xopt,yopt
        integer ifail

        x(1)=x1
        x(2)=x2
        x(3)=x3

        y(1)=y1
        y(2)=y2
        y(3)=y3

        call UTIL_PARABEL(X,Y,A,YP,XOPT,yopt,IFAIL)

        if (ifail.ne.0) then
          print*,'*** Warning in PARABEL: IFAIL .ne. 0'
        endif

      RETURN
      END
