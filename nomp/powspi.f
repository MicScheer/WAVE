*CMZ :  3.00/00 11/03/2013  10.40.59  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.41/13 12/11/2009  16.27.11  by  Michael Scheer
*CMZ :  2.16/08 20/10/2000  11.29.02  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.43.04  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.53.04  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.27  by  Michael Scheer
*-- Author : Michael Scheer
C****************************************************************************
      SUBROUTINE POWSPI(X,Y,IWALL,IMODE,IPOL,ISTART,N)
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

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.

C--- EVALUATES SPLINE-COEFFCIENTS (POWS2). SEE NUMERICAL RECIPIES, PAGE 88

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,spect.
      include 'spect.cmn'
*KEND.


      INTEGER KHI,KLO,K,IPOL,IMODE,IWALL,N,IS,ISTART
      REAL*4 X,Y
      DOUBLE PRECISION XI,XE,HH,AA,BB

      IS=ISTART-1
      XI=RADPOW(IWALL,IPOL,IS+1)
      XE=RADPOW(IWALL,IPOL,IS+N)

      IF (
     &       XE.GT.XI.AND.X.GE.XI.AND.X.LE.XE
     &  .OR.
     &       XI.GT.XE.AND.X.GE.XE.AND.X.LE.XI) THEN

      KLO=IS+1
      KHI=IS+N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(RADPOW(IWALL,IPOL,K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF

      HH=RADPOW(IWALL,IPOL,KHI)-RADPOW(IWALL,IPOL,KLO)
      IF (HH.EQ.0.) THEN
          WRITE(6,*) '*** ERROR IN POWSPI ***'
          STOP
      ENDIF
      AA=(RADPOW(IWALL,IPOL,KHI)-X)/HH
      BB=(X-RADPOW(IWALL,IPOL,KLO))/HH
      Y=   AA*RADPOW(IWALL+2*IMODE+2,IPOL,KLO)
     &    +BB*RADPOW(IWALL+2*IMODE+2,IPOL,KHI)
     &    +((AA**3-AA)*POWS2(KLO-IS)+(BB**3-BB)*POWS2(KHI-IS))
     &    *(HH**2)/6.


      ELSE
          Y=0.
      ENDIF

      RETURN
      END
