*CMZ :  2.53/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  11.03.35  by  Michael Scheer
*CMZ :  2.20/01 12/12/2000  12.09.19  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.12/00 27/05/99  18.01.14  by  Michael Scheer
*CMZ :  2.10/01 07/05/99  12.21.34  by  Michael Scheer
*CMZ : 00.01/10 20/08/96  12.00.15  by  Michael Scheer
*CMZ : 00.01/08 21/06/95  17.17.07  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.22.06  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.48.13  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.17  by  Michael Scheer
*-- Author : Michael Scheer
C*******************************************************************************
      Subroutine MYBMOVE(XI,YI,ZI,VXI,VYI,VZI,BXIN,BYIN,BZIN,DT,
     &             XO,YO,ZO,VXO,VYO,VZO,VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT)
C*******************************************************************************
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
C
C
C     *** BESCHLEUNIGUNG AM ENDE DES INTERVALLS
C
C
C BERECHNET DIE EXAKTE 3-DIM Trajektorie EINES TEILCHENS IN EINEM
C 3-DIM HOMOGENEN Magnetfeld
C
C LAENGEN IN METERN
C GeschwindigkeitEN IN M/SEC
C ZEIT IN SEKUNDEN
C B-FELDER IN TESLA V SEC/M**2
C
C*******************************************************************

      IMPLICIT NONE

*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DOUBLE PRECISION XI,YI,ZI,VXI,VYI,VZI,BX,BY,BZ,DT,
     &                   XO,YO,ZO,VXO,VYO,VZO,VXP,VYP,VZP,GAMMA,
     &                   V,B,EBX,EBY,EBZ,ERX,ERY,ERZ,
     &                   EPX,EPY,EPZ,PPER,RHO,PHI,COSPHI,SINPHI,
     &                   XC,YC,ZC,FRX,FRY,FRZ,FPX,FPY,FPZ,VPAR,VPER,
     &                   V1,B1,OM,APER,ER
C     &                  ,VMX,VMY,VMZ
     &,CONST,BXIN,BYIN,BZIN

        DOUBLE PRECISION BMOVECUT

      INTEGER ICHARGE

        IF (ICHARGE.GT.0) THEN
            BX=-BXIN
            BY=-BYIN
            BZ=-BZIN
        ELSE
            BX=BXIN
            BY=BYIN
            BZ=BZIN
        ENDIF

      V=SQRT(VXI*VXI+VYI*VYI+VZI*VZI)
      B=SQRT(BX*BX+BY*BY+BZ*BZ)

      V1=1.D0/V
      IF (B.GT.BMOVECUT) THEN
         B1=1.D0/B
      ELSE
              XO=XI+VXI*DT
              YO=YI+VYI*DT
              ZO=ZI+VZI*DT
              VXO=VXI
              VYO=VYI
              VZO=VZI
              VXP=0.D0
              VYP=0.D0
              VZP=0.D0
              GOTO 999    !RETURN
      ENDIF

      EBX=BX*B1
      EBY=BY*B1
      EBZ=BZ*B1

      ERX=VYI*EBZ-VZI*EBY
      ERY=VZI*EBX-VXI*EBZ
      ERZ=VXI*EBY-VYI*EBX

      ER=1.D0/SQRT(ERX*ERX+ERY*ERY+ERZ*ERZ)
      ERX=ERX*ER
      ERY=ERY*ER
      ERZ=ERZ*ER

      EPX=EBY*ERZ-EBZ*ERY
      EPY=EBZ*ERX-EBX*ERZ
      EPZ=EBX*ERY-EBY*ERX

C VELOCITY PARALLEL TO B

      VPAR=VXI*EBX+VYI*EBY+VZI*EBZ

C VELOCITY PERPENDICULAR TO B

      VPER=VXI*EPX+VYI*EPY+VZI*EPZ

C MOMENTUM PERPENDICULAR TO B [GEV]

      PPER=GAMMA*EMASSE1*VPER/CLIGHT1
      RHO=PPER/(CLIGHT1*B)
      OM=VPER/RHO
      PHI=OM*DT
      APER=VPER*OM
      COSPHI=COS(PHI)
      SINPHI=SIN(PHI)

      XC=XI-RHO*ERX
      YC=YI-RHO*ERY
      ZC=ZI-RHO*ERZ

      FRX=ERX*COSPHI+EPX*SINPHI
      FRY=ERY*COSPHI+EPY*SINPHI
      FRZ=ERZ*COSPHI+EPZ*SINPHI

      FPX=-ERX*SINPHI+EPX*COSPHI
      FPY=-ERY*SINPHI+EPY*COSPHI
      FPZ=-ERZ*SINPHI+EPZ*COSPHI

      XO=XC+RHO*FRX+VPAR*EBX*DT
      YO=YC+RHO*FRY+VPAR*EBY*DT
      ZO=ZC+RHO*FRZ+VPAR*EBZ*DT

      VXO=VPER*FPX+VPAR*EBX
      VYO=VPER*FPY+VPAR*EBY
      VZO=VPER*FPZ+VPAR*EBZ

C     VMX=(VXI+VXO)*0.5D0
C     VMY=(VYI+VYO)*0.5D0
C     VMZ=(VZI+VZO)*0.5D0

C     CONST=-ECHARGE1/(GAMMA*EMASSKG1)
C     VXP=(VMY*BZ-VMZ*BY)*CONST
C     VYP=(VMZ*BX-VMX*BZ)*CONST
C     VZP=(VMX*BY-VMY*BX)*CONST

      CONST=-ECHARGE1/(GAMMA*EMASSKG1)
      VXP=(VYO*BZ-VZO*BY)*CONST
      VYP=(VZO*BX-VXO*BZ)*CONST
      VZP=(VXO*BY-VYO*BX)*CONST


999   RETURN
      END
