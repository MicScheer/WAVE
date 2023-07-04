*CMZ :  4.00/11 17/05/2021  11.31.37  by  Michael Scheer
*CMZ :  3.01/03 19/03/2014  12.24.14  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/19 07/06/2011  14.08.31  by  Michael Scheer
*CMZ :  2.63/03 07/05/2008  14.17.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.14.32  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  18.41.04  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.55.02  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.59  by  Michael Scheer
*-- Author : Michael Scheer
C************************************************************************
      SUBROUTINE SPLINTBz(N,X,Y)
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

      INTEGER NOLD,N,KLO,KHI,KLOLD,K
      DOUBLE PRECISION Y,X,XA1OLD,XANOLD,H,A,B

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEND.
      DOUBLE PRECISION XA(NBTABP),YA(NBTABP),Y2A(NBTABP)

      COMMON/BTABCz/XA,YA,Y2A

      save klold,nold,xa1old,xanold

      DATA KLOLD/1/,NOLD/-99/
      DATA XA1OLD/-9999.D0/,XANOLD/-9999./

      IF(     XA(1).LT.XA(N).AND.(X.LT.XA(1).OR.X.GT.XA(N))
     &  .OR.
     &  XA(N).LT.XA(1).AND.(X.LT.XA(N).OR.X.GT.XA(1)))
     &  STOP '***SR SPLINTBz: X OUT OF RANGE ***'

      IF (NOLD.ne.N) then
        klo=1
      else IF (
     &    XA(1).EQ.XA1OLD
     &    .AND. XA(N).EQ.XANOLD
     &    .AND. X.GT.XA(KLOLD)
     &    ) THEN
        KLO=KLOLD
      ELSE
        KLO=1
      ENDIF

      IF (X.LT.XA(KLO+1)) THEN
      KHI=KLO+1
      GOTO 2
      ENDIF

      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF

2     H=XA(KHI)-XA(KLO)
      IF (H.LE.0.) THEN
        WRITE(6,*) 'Bad XA input.'
        STOP
      ENDIF
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.

      KLOLD=KLO
      NOLD=N
      XA1OLD=XA(1)
      XANOLD=XA(N)


      RETURN
      END
