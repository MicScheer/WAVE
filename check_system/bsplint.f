*CMZ :  3.02/00 24/09/2014  13.51.08  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.66/19 07/06/2011  14.08.31  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  16.46.16  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.48.34  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.15  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BSPLINT(XA,YAX,YAY,Y2AX,Y2AY,N,X,YX,YY,KOLD,KLO)
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

      INTEGER N,INC,KOLD,KLO,KHI,K
      REAL XA,YAX,YAY,Y2AX,Y2AY,X,YX,YY,B
      REAL H,A

      DIMENSION XA(*),YAX(*),YAY(*),Y2AX(*),Y2AY(*)

C SUBROUTINE MACHT SPLINE-INTERPOLATIONS. SIEHE "NUMERICAL RECIPIES"

      IF(X.GT.XA(N).OR.X.LT.XA(1)) THEN
          STOP '*** S/R BSPLINT: X OUT OF RANGE ***'
      ELSE IF (X.EQ.XA(N)) THEN
          YX=YAX(N)
          YY=YAY(N)
          RETURN
      ELSE IF (X.EQ.XA(1)) THEN
          YX=YAX(1)
          YY=YAY(1)
          RETURN
      ENDIF

      KLO=KOLD

      IF(KLO.LT.1 .OR. KLO.GT.N) THEN
          KLO=1
          KHI=N
          GOTO 3
      ENDIF

      INC=1
      IF (X.GT.XA(KLO)) THEN
1         KHI=KLO+INC
          IF (KHI.GT.N) THEN
         KHI=N
          ELSE IF (X.GE.XA(KHI)) THEN
         KLO=KHI
         INC=INC+INC
         GOTO 1
          ENDIF
      ELSE
          KHI=KLO
2         KLO=KHI-INC
          IF (KLO.LT.1) THEN
         KLO=1
          ELSE IF (X.LT.XA(KLO)) THEN
         KHI=KLO
         INC=INC+INC
         GOTO 2
          ENDIF
      ENDIF
3     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 3
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.LE.0.    .OR. X.LT.XA(KLO).OR.X.GT.XA(KHI)) THEN
        WRITE(6,*) '*** SR BSPLINT: Bad XA input. ***'
        STOP
      ENDIF
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      YX=A*YAX(KLO)+B*YAX(KHI)+
     *      ((A**3-A)*Y2AX(KLO)+(B**3-B)*Y2AX(KHI))*(H**2)/6.
      YY=A*YAY(KLO)+B*YAY(KHI)+
     *      ((A**3-A)*Y2AY(KLO)+(B**3-B)*Y2AY(KHI))*(H**2)/6.
      RETURN
      END
