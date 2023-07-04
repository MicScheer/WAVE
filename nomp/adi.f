*CMZ :  2.70/12 01/03/2013  15.45.11  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.02.51  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.58.51  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  14.02.13  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.46.30  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.22  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE ADI(GAMMA,B0,XLAM0,FASYM,F0,DI2AHW,DI5AHW,BETAHW,
     &                 B0PK,XLPK,BETADIK,FB0M,F0P,FBETP,
     &                 DI2ADIK,DI5ADIK,CHI2MIN)

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

C--- BERECHENET FUER DEN ASYM. HALBACH-WLS EIN ASYM. DIPOLAEQUIVALENT
C    MIT GLEICHER ASYMMETRIE FB0M=FB0M(FASYM)
C    FORMFAKTOREN WERDEN AUS NUMERISCHEN BERECHNETEN TABELLEN INTERPOLIERT
      IMPLICIT NONE

      INTEGER IX,IB,NB0P,NXLP
      DOUBLE PRECISION B0,XLAM0,FASYM,F0,B0P,XLP,FB0M,FB0MFUN,FBETP,F0P,FM,F0PFUN,
     &         FBETPFUN
      DOUBLE PRECISION E,GAMMA,CHI2,CHI2MIN
      DOUBLE PRECISION B0PMX,B0PMN,DB0P,RHOP,B0PK,PHIP
      DOUBLE PRECISION XLPMX,XLPMN,DXLP,XLPK
      DOUBLE PRECISION FN,BETADI,BETADIK,BETAHW
      DOUBLE PRECISION DI2AHW,DI5AHW,DI2ADI,DI5ADI,DI2ADIK,DI5ADIK

C     DATA NB0P/100/
C     DATA NXLP/100/

      F0=F0

      E=GAMMA*511003.3732832001D-9

      B0PMX=2.*B0
      B0PMN=0.
      DB0P=0.05

      XLPMX=3.*XLAM0
      XLPMN=0.
      DXLP=0.01

      FB0M=FB0MFUN(FASYM) !EMPIRISCHER FIT BZW INTERPOLATION
      FM=DSQRT(1+1./FB0M**2)**3
      FN=DSQRT(((1.+FASYM)**2/(8.*FASYM)))**3
      F0P=F0PFUN(FASYM)
      FBETP=FBETPFUN(FASYM)

      NB0P=NINT((B0PMX-B0PMN)/DB0P)
      NXLP=NINT((B0PMX-B0PMN)/DXLP)

      CHI2MIN=1.D30
      DO IB=1,NB0P
      DO IX=1,NXLP

          B0P=B0PMN+DB0P*IB
          RHOP=E/B0P/0.3
          XLP=XLPMN+DXLP*IX
          PHIP=XLP/RHOP

          DI2ADI=XLP/(4.*RHOP**2)*(1.+1./FB0M**2)*(1.+1./FASYM)
          DI5ADI=DI2ADI*F0P*FM*FN*PHIP**3
          BETADI=FBETP*XLP*(1.+FASYM)/2.

C         CHI2=((BETAHW-BETADI)/BETAHW)**2+
C     &           ((DI2AHW-DI2ADI)/DI2AHW)**2+
C     &      ((DI5AHW-DI5ADI)/DI5AHW)**2

          CHI2=((BETAHW-BETADI)/BETAHW)**2+
     &            ((DI2AHW-DI2ADI)/DI2AHW)**2+
     &       ((DI5AHW/DI2AHW-DI5ADI/DI2ADI)/(DI5AHW/DI2AHW))**2


          IF (CHI2.LT.CHI2MIN) THEN

         B0PK=B0P
         XLPK=XLP
         DI2ADIK=DI2ADI
         DI5ADIK=DI5ADI
         BETADIK=BETADI
         CHI2MIN=CHI2

          ENDIF


      ENDDO
      ENDDO

      IF(B0PK.EQ.B0PMN.OR.B0PK.EQ.B0PMX)
     &     STOP'*** S/R ADI: B0P-BOUNDARY TOUCHED'
      IF(XLPK.EQ.XLPMN.OR.XLPK.EQ.XLPMX)
     &     STOP'*** S/R ADI: XLP-BOUNDARY TOUCHED'

      CALL ADIBETA(GAMMA,B0PK,XLPK,FASYM,FB0M)

      RETURN
      END
