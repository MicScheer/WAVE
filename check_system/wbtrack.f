*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.63/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.61/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  10.48.09  by  Michael Scheer
*CMZ :  2.47/12 24/06/2003  15.45.35  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  11.29.39  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  2.12/00 02/06/99  13.56.19  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.14.52  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.14  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WBTRACK(X1,Y1,Z1,VX1,VY1,VZ1,
     &  XF0,YF0,ZF0,VXF0,VYF0,VZF0,
     &  X2,Y2,Z2,VX2,VY2,VZ2,DTIM0,GAMMA,IERROR)
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
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEND.

      INTEGER IZAEHL,IERROR

      DOUBLE PRECISION X1,Y1,Z1,VX1,VY1,VZ1,
     &  XF0,YF0,ZF0,VXF0,VYF0,VZF0,
     &  X2,Y2,Z2,VX2,VY2,VZ2,DTIM0,GAMMA,DSMAX
      DOUBLE PRECISION EWSFX,EWSFY,EWSFZ,V0,VXPINT,VYPINT,VZPINT
      DOUBLE PRECISION BX2,BY2,BZ2,X2B,Y2B,Z2B
      DOUBLE PRECISION AX2,AY2,AZ2,VXP,VYP,VZP
      DOUBLE PRECISION DIST2,DIST1,DT,DDT
     &  ,DGAMMA

      IERROR=0
      V0=DSQRT(VX1**2+VY1**2+VZ1**2)
      EWSFX=VXF0/V0
      EWSFY=VYF0/V0
      EWSFZ=VZF0/V0

C030293  DSMAX=DTIM0*CLIGHT1*1.D-8
      DSMAX=DMAX1(DTIM0*CLIGHT1*1.D-8,1.D-10)

      DT=DTIM0/2.D0

      X2=X1
      Y2=Y1
      Z2=Z1

      VX2=VX1
      VY2=VY1
      VZ2=VZ1

         IZAEHL=0

1000     IZAEHL=IZAEHL+1

         X1=X2
         Y1=Y2
         Z1=Z2

         VX1=VX2
         VY1=VY2
         VZ1=VZ2

         X2B=X1+VX1*DT
         Y2B=Y1+VY1*DT
         Z2B=Z1+VZ1*DT

         CALL MYBFELD(X2B,Y2B,Z2B,BX2,BY2,BZ2,AX2,AY2,AZ2)

         CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DTIM0,
     &            X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

C EWSF IS NORMAL VECTOR OF PERPENDICULARE PLANE AT THE END OF THE REFERENCE ORBI
C DIST IS DISTANCE OF ELECTRON TO THIS PLANE
C TRACKING STOPS IF TRAJECTORIE HITS THIS PLANE

         DIST2=(X2-XF0)*EWSFX+(Y2-YF0)*EWSFY+(Z2-ZF0)*EWSFZ

      IF ( DIST2 .LT. 0.0  )  GOTO 1000


C--- ENDE OF TRAJECTORY, DIST2 NOT EXACTLY ZERO, CORRECT X2

      DIST1=(X1-XF0)*EWSFX+(Y1-YF0)*EWSFY+(Z1-ZF0)*EWSFZ
      DDT=DTIM0*DABS(DIST1)/(DABS(DIST1)+DABS(DIST2))

      CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDT,
     &  X2,Y2,Z2,VX2,VY2,VZ2,
     &  VXPINT,VYPINT,VZPINT,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

      DIST2=(X2-XF0)*EWSFX+(Y2-YF0)*EWSFY+(Z2-ZF0)*EWSFZ

C030293  IF (DIST2.GT.DSMAX) THEN
      IF (DABS(DIST2).GT.DSMAX) THEN
          IERROR=1
      ENDIF

      END
