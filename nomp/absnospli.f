*CMZ :  3.02/00 24/09/2014  13.51.08  by  Michael Scheer
*CMZ :  3.01/03 19/03/2014  12.24.14  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  15.45.11  by  Michael Scheer
*CMZ :  2.66/19 07/06/2011  14.08.31  by  Michael Scheer
*CMZ :  2.36/00 07/11/2001  14.33.03  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.21.38  by  Michael Scheer
*CMZ : 00.01/01 21/06/94  20.13.12  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.46.24  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.34  by  Michael Scheer
*-- Author : Michael Scheer

      SUBROUTINE ABSNOSPLI(XA,YA,N,X,Y,IERR,IMODE)

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
      INTEGER NOLD,N,KLO,KHI,KLOLD,K,IERR,IMODE
      DOUBLE PRECISION Y,X,XA1OLD,XANOLD,H,YA1OLD,YANOLD

      DOUBLE PRECISION XA(*),YA(*),BB

      save klold,nold,xa1old,xanold,ya1old,yanold

      DATA KLOLD/1/,NOLD/-99/
      DATA XA1OLD/-9999.D0/,XANOLD/-9999./
      DATA YA1OLD/-9999.D0/,YANOLD/-9999./

      IF( XA(1).GE.XA(N)) THEN
        WRITE(6,*)
        WRITE(6,*) '***ERROR IN ABSNOSPLI: X NOT INCREASING ***'
        WRITE(6,*)
        IERR=1
        RETURN
      ENDIF

      IF( XA(1).LT.XA(N).AND.(X.LT.XA(1).OR.X.GT.XA(N))
     &    .OR.
     &    XA(N).LT.XA(1).AND.(X.LT.XA(N).OR.X.GT.XA(1))) THEN
        WRITE(6,*)
        WRITE(6,*) '***ERROR IN ABSNOSPLI: X OUT OF RANGE ***'
        WRITE(6,*)
        WRITE(6,*) 'X-RANGE:',XA(1),'-',XA(N)
        WRITE(6,*) 'X:      ',X
        WRITE(6,*)
        IERR=1
        RETURN
      ENDIF

      IF (NOLD.EQ.N) THEN
        IF (
     &      XA(1).EQ.XA1OLD
     &     .AND. XA(N).EQ.XANOLD
     &     .AND. YA(1).EQ.YA1OLD
     &     .AND. YA(N).EQ.YANOLD
     &     .AND. X.GT.XA(KLOLD)
     &   ) THEN
           KLO=KLOLD
        ELSE
           KLO=1
        ENDIF
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

      IF (IMODE.GT.0) THEN

C INTERPOLATION OF Y BY Y=AA*X**BB

        IF (H.NE.0.) THEN
          BB=DLOG(YA(KHI)/YA(KLO))/DLOG(XA(KHI)/XA(KLO))
          Y=YA(KLO)*DEXP(BB*(DLOG(X/XA(KLO))))
        ELSE
C           IF (YA(KLO).NE.YA(KHI)) STOP '*** SR ABSNOSPLI: Bad Input ***'
          Y=YA(KLO)
        ENDIF

      ELSE !(IMODE)

C LINEAR INTERPOLATION

        IF (H.NE.0.) THEN
          Y=YA(KLO)+(X-XA(KLO))/H*(YA(KHI)-YA(KLO))
        ELSE
          Y=YA(KLO)
        ENDIF

      ENDIF !(IMODE)

      KLOLD=KLO
      NOLD=N
      XA1OLD=XA(1)
      XANOLD=XA(N)
      YA1OLD=YA(1)
      YANOLD=YA(N)

      RETURN
      END
