*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.63/05 22/07/2009  07.43.29  by  Michael Scheer
*CMZ :  2.48/04 12/03/2004  15.40.31  by  Michael Scheer
*CMZ :  2.47/07 14/04/2003  15.17.05  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  16.33.46  by  Michael Scheer
*CMZ : 00.01/04 28/11/94  18.36.50  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.55  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE FSPLINDX(DX,Y,N,YP1,YPN,Y2)
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
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEND.

      INTEGER N,J
      DOUBLE PRECISION DX,Y(N),Y2(N)

      DOUBLE PRECISION  C(NDOBSVZP+NDOBSVYP)
      DOUBLE PRECISION AA(NDOBSVZP+NDOBSVYP)
      DOUBLE PRECISION BB(NDOBSVZP+NDOBSVYP)
      DOUBLE PRECISION CC(NDOBSVZP+NDOBSVYP)


      DOUBLE PRECISION YP1,YPN,H6,H3

      Y2(1)=YP1
      Y2(N)=YPN

      IF (N.LT.3) RETURN

      C(1)=YP1
      C(N)=YPN

      H6=1./6.D0
      H3=4.D0*H6

      BB(1)=1.D0
      CC(1)=0.D0


      DO J=2,N-1
          AA(J)=H6
          BB(J)=H3
          CC(J)=H6
          C(J)=(Y(J+1)-2.D0*Y(J)+Y(J-1))/DX**2
      ENDDO !J

      DO J=2,N-1

          BB(J)=BB(J)-AA(J)*CC(J-1)
           C(J)= C(J)-AA(J)* C(J-1)
C030414          AA(J)=AA(J)-AA(J)*BB(J-1)

          CC(J)=CC(J)/BB(J)
           C(J)= C(J)/BB(J)
          BB(J)=1.D0

      ENDDO !J

      DO J=N-1,2,-1
         Y2(J)=C(J)-CC(J)*Y2(J+1)
      ENDDO

      RETURN
      END
