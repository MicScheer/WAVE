*CMZ :  4.00/15 05/04/2022  11.56.08  by  Michael Scheer
*CMZ :  3.04/00 04/01/2018  18.10.43  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/20 06/07/2011  09.41.56  by  Michael Scheer
*CMZ :  2.63/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.61/02 14/03/2007  15.48.38  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  10.55.32  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  11.25.52  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.12/00 02/06/99  13.56.19  by  Michael Scheer
*CMZ :  1.02/00 19/12/97  17.09.57  by  Michael Scheer
*CMZ : 00.02/01 18/12/96  11.40.20  by  Michael Scheer
*CMZ : 00.01/10 16/07/96  14.53.50  by  Michael Scheer
*CMZ : 00.01/09 11/04/96  17.41.50  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  09.51.12  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.41.55  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.31  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE TRACKBMAG(ISNORDER,X1,Y1,Z1,VX1,VY1,VZ1,
     &                      XF0,YF0,ZF0,EWSFX,EWSFY,EWSFZ,
     &  X2,Y2,Z2,VX2,VY2,VZ2,DTIM,BSHIFT,GAMMA,IMag,
     &  BMOVECUT,iustep,ieneloss)
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

C--- CLONE OF TRACKSHORT, BUT CALLING BDI or BDH INSTEAD OF MYBFELD


C---     INPUT:
C     ISNORDER   order of tracking
C     X1,Y1,Z1    coordinates of electron at the starting point
C     VX1,VY1,VZ1    velocity at the starting point
C     XF0,YF0,ZF0    point of the plane where tracking stops
C     EWSFX,EWSFY,EWSFZ normal vector of the plane
C     DTIM        time intervall of one tracking step
C     BSHIFT         fraction of step where magnetic field is
C              determined for this step
C     GAMMA       relativistic factor

C---  OUTPUT:
C     X2,Y2,Z2    final coordinates of the electron
C     VX2,VY2,VZ2    final velocity of the electron

         IMPLICIT NONE

      INTEGER ISNORDER
      INTEGER IM,IUSTEP,IENELOSS,imag

         DOUBLE PRECISION X1,Y1,Z1,VX1,VY1,VZ1,X2,Y2,Z2,VX2,VY2,VZ2
     &           ,DTIM,BSHIFT,X2B,Y2B,Z2B,BX1,BY1,BZ1,BX2,BY2,BZ2
     &           ,GAMMA,DT,VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM
     &           ,X2INT,Y2INT,Z2INT,DDT,DDDT,DDT2
     &           ,VX2INT,VY2INT,VZ2INT,VXPINT,VYPINT,VZPINT
     &           ,VXP,VYP,VZP,XOLD,YOLD,ZOLD,VXOLD,VYOLD,VZOLD
     &           ,X3INT,Y3INT,Z3INT,VX3INT,VY3INT,VZ3INT,DDDDT,DDDDT2
     &           ,EWSFX,EWSFY,EWSFZ,XF0,YF0,ZF0,DIST1,DIST2,DISTI,BMOVECUT

      DOUBLE PRECISION T1,T2
     &  ,DGAMMA

*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      IF (IENELOSS.NE.0) THEN
        WRITE(16,*)
        WRITE(16,*)
     &    ' *** WARNING IN TRACKBMAG: IENELOSS .NE. 0, NOT YET IMPLEMENTED ***'
        WRITE(16,*)
     &    ' *** BE CAREFUL!! ***'
        WRITE(16,*)
        WRITE(6,*)
        WRITE(6,*)
     &    ' *** WARNING IN TRACKBMAG: IENELOSS .NE. 0, NOT YET IMPLEMENTED ***'
        WRITE(6,*)
     &    ' *** BE CAREFUL!! ***'
        WRITE(6,*)
      ENDIF

      im=iabs(imag)

      DT=DTIM*BSHIFT

      T1=0.0

      BX1=0.D0
      BY1=0.D0
      BZ1=0.D0

      XOLD=X1
      YOLD=Y1
      ZOLD=Z1
      VXOLD=VX1
      VYOLD=VY1
      VZOLD=VZ1

      X2=X1
      Y2=Y1
      Z2=Z1
      T2=T1

      VX2=VX1
      VY2=VY1
      VZ2=VZ1

      BX2=BX1
      BY2=BY1
      BZ2=BZ1

C--- LOOP DER TRAJEKTORIE

1000  X1=X2
      Y1=Y2
      Z1=Z2

      T1=T2

      VX1=VX2
      VY1=VY2
      VZ1=VZ2

      BX1=BX2
      BY1=BY2
      BZ1=BZ2

      IF (ISNORDER.EQ.0) THEN

        X2B=X1+VX1*DT
        Y2B=Y1+VY1*DT
        Z2B=Z1+VZ1*DT

      ELSE

        CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX1,BY1,BZ1,DT,X2B,Y2B,Z2B,
     &    VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM,GAMMA,ICHARGE,BMOVECUT,
     &    IUSTEP,IENELOSS,DGAMMA)

      ENDIF

      if (imag.gt.0) then
        CALL BDI(X2B,Y2B,Z2B,BX2,BY2,BZ2,IM)
      else
        CALL BDH(X2B,Y2B,Z2B,BX2,BY2,BZ2,IM)
      endif

      CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DTIM,
     &  X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,
     &  IUSTEP,IENELOSS,DGAMMA)

      T2=T1+DTIM

C EWSF IS NORMAL VECTOR OF PERPENDICULARE PLANE AT THE END OF THE REFERENCE ORBI
C DIST IS DISTANCE OF ELECTRON TO THIS PLANE
C TRACKING STOPS IF TRAJECTORIE HITS THIS PLANE

      DIST2=(X2-XF0)*EWSFX+(Y2-YF0)*EWSFY+(Z2-ZF0)*EWSFZ

      IF ( DIST2 .LT. 0.0  )  GOTO 1000

C--- ENDE OF TRAJECTORY, DIST2 NOT EXACTLY ZERO, CORRECT X2

      DIST1=(X1-XF0)*EWSFX+(Y1-YF0)*EWSFY+(Z1-ZF0)*EWSFZ

      DDT=DTIM*DABS(DIST1)/(DABS(DIST1)+DABS(DIST2))
      DDT2=DDT*BSHIFT

      IF(ISNORDER.EQ.0) THEN

         X2B=X1+VX1*DDT2
         Y2B=Y1+VY1*DDT2
         Z2B=Z1+VZ1*DDT2

      ELSE

         CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX1,BY1,BZ1,DDT2,X2B,Y2B,Z2B,
     &     VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

       ENDIF

       if (imag.gt.0) then
         CALL BDI(X2B,Y2B,Z2B,BX2,BY2,BZ2,IM)
       else
         CALL BDH(X2B,Y2B,Z2B,BX2,BY2,BZ2,IM)
       endif

       CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDT,
     &   X2INT,Y2INT,Z2INT,VX2INT,VY2INT,VZ2INT,
     &   VXPINT,VYPINT,VZPINT,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

       DISTI=(X2INT-XF0)*EWSFX+(Y2INT-YF0)*EWSFY+(Z2INT-ZF0)*EWSFZ
       DDDT=DDT
       IF (DIST1.NE.0.) DDDT=DDT-DDT*DISTI/DABS(DIST1)

       CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDDT,
     &   X3INT,Y3INT,Z3INT,VX3INT,VY3INT,VZ3INT,
     &   VXPINT,VYPINT,VZPINT,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

       DISTI=(X3INT-XF0)*EWSFX+(Y3INT-YF0)*EWSFY+(Z3INT-ZF0)*EWSFZ
       DDDDT=DDDT
       IF (DIST1.NE.0.) DDDDT=DDDT-DDDT*DISTI/DABS(DIST1)
       DDDDT2=BSHIFT*DDDDT

       CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDDDT,
     &   X2,Y2,Z2,VX2,VY2,VZ2,
     &   VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

       T2=T1+DDDDT

       X1=XOLD
       Y1=YOLD
       Z1=ZOLD
       VX1=VXOLD
       VY1=VYOLD
       VZ1=VZOLD

      RETURN
      END
