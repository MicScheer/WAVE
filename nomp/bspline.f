*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.13/09 08/03/2000  18.08.37  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.14  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BSPLINE(X,YX,YY,N,YP1,YPN,Y2X,Y2Y,UX,UY)
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

      INTEGER K,I,N

      DOUBLE PRECISION X,YX,YY,Y2X,Y2Y,UX,UY,QN,UNX,UNY,SIG,PX,PY,YPN,YP1
        DIMENSION X(N),YX(N),YY(N),Y2X(N),Y2Y(N),UX(N),UY(N)

      IF (YP1.GT..99E30) THEN
        Y2X(1)=0.
        UX(1)=0.
        Y2Y(1)=0.
        UY(1)=0.
      ELSE
        Y2X(1)=-0.5
        UX(1)=(3./(X(2)-X(1)))*((YX(2)-YX(1))/(X(2)-X(1))-YP1)
        Y2Y(1)=-0.5
        UY(1)=(3./(X(2)-X(1)))*((YY(2)-YY(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        PX=SIG*Y2X(I-1)+2.
        Y2X(I)=(SIG-1.)/PX
        UX(I)=(6.*((YX(I+1)-YX(I))/(X(I+1)-X(I))-(YX(I)-YX(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*UX(I-1))/PX
        PY=SIG*Y2Y(I-1)+2.
        Y2Y(I)=(SIG-1.)/PY
        UY(I)=(6.*((YY(I+1)-YY(I))/(X(I+1)-X(I))-(YY(I)-YY(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*UY(I-1))/PY
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UNX=0.
        UNY=0.
      ELSE
        QN=0.5
        UNX=(3./(X(N)-X(N-1)))*(YPN-(YX(N)-YX(N-1))/(X(N)-X(N-1)))
        UNY=(3./(X(N)-X(N-1)))*(YPN-(YY(N)-YY(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2X(N)=(UNX-QN*UX(N-1))/(QN*Y2X(N-1)+1.)
      Y2Y(N)=(UNY-QN*UY(N-1))/(QN*Y2Y(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2X(K)=Y2X(K)*Y2X(K+1)+UX(K)
        Y2Y(K)=Y2Y(K)*Y2Y(K+1)+UY(K)
12    CONTINUE
      RETURN
      END
