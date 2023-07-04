*CMZ :  4.00/16 22/07/2022  09.45.26  by  Michael Scheer
*CMZ :  4.00/15 29/04/2022  08.25.50  by  Michael Scheer
*CMZ :  3.07/00 05/03/2019  13.35.52  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.67/02 16/03/2012  16.36.45  by  Michael Scheer
*CMZ :  2.66/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.65/02 28/09/2009  12.44.41  by  Michael Scheer
*CMZ :  2.63/05 10/08/2009  20.48.43  by  Michael Scheer
*CMZ :  2.61/02 14/03/2007  15.48.38  by  Michael Scheer
*CMZ :  2.53/04 09/02/2005  11.06.23  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  10.54.01  by  Michael Scheer
*CMZ :  2.47/17 12/09/2003  10.04.22  by  Michael Scheer
*CMZ :  2.47/12 24/06/2003  16.20.58  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  11.26.52  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.12/00 02/06/99  13.56.19  by  Michael Scheer
*CMZ : 00.02/01 18/12/96  11.40.20  by  Michael Scheer
*CMZ : 00.01/10 16/07/96  14.53.50  by  Michael Scheer
*CMZ : 00.01/09 11/04/96  17.41.50  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  09.51.12  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.41.55  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.31  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE TRACKSHORT(ISNORDER,Xi,Yi,Zi,VXi,VYi,VZi,
     &  XF0,YF0,ZF0,EWSFX,EWSFY,EWSFZ,
     &  X2,Y2,Z2,t2,VX2,VY2,VZ2,DTIM,BSHIFT,GAMMA0,BMOVECUT
     &  ,IUSTEP,IENELOSS,GAMMAL)
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

C---     INPUT:
C     ISNORDER   order of tracking
C     XI,YI,ZI    coordinates of electron at the starting point
C     VXI,VYI,VZI    velocity at the starting point
C     XF0,YF0,ZF0    point of the plane where tracking stops
C     EWSFX,EWSFY,EWSFZ normal vector of the plane
C     DTIM        time intervall of one tracking step
C     BSHIFT         fraction of step where magnetic field is
C              determined for this step
C     GAMMA       relativistic factor

C---  OUTPUT:
C     X2,Y2,Z2,T2    final coordinates and time-of-flight of the electron
C     VX2,VY2,VZ2    final velocity of the electron

      IMPLICIT NONE

      INTEGER ISNORDER,I,IROI1,IROI2,IUSTEP,IENELOSS,IWARN

      DOUBLE PRECISION X1,Y1,Z1,VX1,VY1,VZ1,X2,Y2,Z2,VX2,VY2,VZ2
     &  ,DTIM,BSHIFT,X2B,Y2B,Z2B,BX1,BY1,BZ1,BX2,BY2,BZ2
     &  ,DGAMSUM,GAMMA,GAMMAL,GAMMA0,DT,VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM
     &  ,X2INT,Y2INT,Z2INT,DDT,DDDT,DDT2
     &  ,VX2INT,VY2INT,VZ2INT,VXPINT,VYPINT,VZPINT
     &  ,VXP,VYP,VZP,XOLD,YOLD,ZOLD,VXOLD,VYOLD,VZOLD
     &  ,X3INT,Y3INT,Z3INT,VX3INT,VY3INT,VZ3INT,DDDDT,DDDDT2
     &  ,EWSFX,EWSFY,EWSFZ,XF0,YF0,ZF0,DIST1,DIST2,DISTI
     &  ,AX1,AY1,AZ1,AX2,AY2,AZ2,vxi,vyi,vzi,xi,yi,zi
     &  ,X2BOUND,Y2BOUND,Z2BOUND,VX2BOUND,VY2BOUND,VZ2BOUND
     &  ,X1SAV,Y1SAV,Z1SAV,VX1SAV,VY1SAV,VZ1SAV,BX1SAV,BY1SAV,BZ1SAV
     &  ,AX1SAV,AY1SAV,AZ1SAV,T1SAV,X2SAV,DTIM0,BMOVECUT,BETA,VN

      DOUBLE PRECISION T1,T2,DGAMMA,VXSIGN

*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEND.

      integer :: idebug=0,ical=0

      DATA IWARN/0/

C VXSIGN takes care for the direction of flight, since particle must gain
c energy if tracked back

      x1=xi
      y1=yi
      z1=zi

      vx1=vxi
      vy1=vyi
      vz1=vzi

      IF (VX1.LT.0) THEN
        VXSIGN=-1.0D0
      ELSE
        VXSIGN=1.0D0
      ENDIF

      GAMMA=GAMMA0
      DGAMSUM=0.0D0

      VN=SQRT(VX1*VX1+VY1*VY1+VZ1*VZ1)

      DT=DTIM*BSHIFT
      DTIM0=DTIM

      T1=0.0

      BX1=0.D0
      BY1=0.D0
      BZ1=0.D0

      AX1=0.D0
      AY1=0.D0
      AZ1=0.D0

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

      IROI1=0
      IROI2=1

      DO I=1,NROI
        IF(X1.GE.ROIX(I)) THEN
          IROI1=I
          IROI2=IROI1+1
          GOTO 1000
        ENDIF
      ENDDO

C--- LOOP DER TRAJEKTORIE

1000  X1=X2
      Y1=Y2
      Z1=Z2

      T1=T2

      VX1=VX2
      VY1=VY2
      VZ1=VZ2

      AX1=AX2
      AY1=AY2
      AZ1=AZ2

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

      CALL MYBFELD(X2B,Y2B,Z2B,BX2,BY2,BZ2,AX2,AY2,AZ2)

      CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DTIM,
     &  X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,
     &  IUSTEP,IENELOSS,DGAMMA)

      if (idebug.ne.0) then
        ical=ical+1
        write(77,*) ical,x1,z1,by2
      endif

C--- BOUNDARY CROSSING {

      IF (X2.GT.ROIX(IROI2)) THEN

        IF (IWARN.EQ.0.AND.DTIM.LT.0.0D0) THEN
          WRITE(16,*)
          WRITE(16,*)' *** WARNING IN TRACKSHORT:'
          WRITE(16,*)
     &      ' *** DTIM negative, boundary is probably treated wrong! ***'
          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*)' *** WARNING IN TRACKSHORT:'
          WRITE(6,*)
     &      ' *** DTIM negative, boundary is probably treated wrong! ***'
          WRITE(6,*)
          IWARN=1
        ENDIF

        X1SAV=X1
        Y1SAV=Y1
        Z1SAV=Z1

        T1SAV=T1

        VX1SAV=VX1
        VY1SAV=VY1
        VZ1SAV=VZ1

        AX1SAV=AX1
        AY1SAV=AY1
        AZ1SAV=AZ1

        BX1SAV=BX1
        BY1SAV=BY1
        BZ1SAV=BZ1

        X2SAV=X2

        DIST1=X1-ROIX(IROI2)
        DIST2=X2-ROIX(IROI2)

        DDT=DTIM*DABS(DIST1)/(DABS(DIST1)+DABS(DIST2))
        DDT2=DDT*BSHIFT

        IF(ISNORDER.EQ.0) THEN

          X2B=X1+VX1*DDT2
          Y2B=Y1+VY1*DDT2
          Z2B=Z1+VZ1*DDT2

        ELSE

          CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX1,BY1,BZ1,DDT2,
     &      X2B,Y2B,Z2B,
     &      VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM,GAMMA,ICHARGE,BMOVECUT,
     &      IUSTEP,IENELOSS,DGAMMA)

        ENDIF

        CALL MYBFELD(X2B,Y2B,Z2B,BX2,BY2,BZ2,AX2,AY2,AZ2)

        CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDT,
     &    X2INT,Y2INT,Z2INT,VX2INT,VY2INT,VZ2INT,
     &    VXPINT,VYPINT,VZPINT,GAMMA,ICHARGE,BMOVECUT,
     &    IUSTEP,IENELOSS,DGAMMA)

        DISTI=X2INT-ROIX(IROI2)

        DDDT=DDT
        IF (DIST1.NE.0.) DDDT=DDT-DDT*DISTI/DABS(DIST1)

        CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDDT,
     &    X3INT,Y3INT,Z3INT,VX3INT,VY3INT,VZ3INT,
     &    VXPINT,VYPINT,VZPINT,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

        DISTI=X3INT-ROIX(IROI2)
        DDDDT=DDDT
        IF (DIST1.NE.0.) DDDDT=DDDT-DDDT*DISTI/DABS(DIST1)
        DDDDT2=BSHIFT*DDDDT

        CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDDDT,
     &    X2BOUND,Y2BOUND,Z2BOUND,VX2BOUND,VY2BOUND,VZ2BOUND,
     &    VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

        IF (IENELOSS.NE.0) THEN
          DGAMSUM=DGAMSUM+VXSIGN*DGAMMA
          IF (ABS(DGAMSUM).GT.GAMMA*1.0D-8) THEN
            GAMMA=GAMMA+DGAMSUM
            DGAMSUM=0.0D0
          ENDIF
          BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
          VN=SQRT(VX2BOUND*VX2BOUND+VY2BOUND*VY2BOUND+VZ2BOUND*VZ2BOUND)
          VX2BOUND=VX2BOUND/VN*CLIGHT1*BETA
          VY2BOUND=VY2BOUND/VN*CLIGHT1*BETA
          VZ2BOUND=VZ2BOUND/VN*CLIGHT1*BETA
        ENDIF

C NOW WE ARE ON THE BOUNDARY (X2BOUND) AND PRECEED TO OLD X2 (X2SAV)

        DTIM=DTIM*(X2SAV-X2BOUND)/(X2SAV-X1SAV)
        DT=DTIM*BSHIFT

        X1=X2BOUND
        Y1=Y2BOUND
        Z1=Z2BOUND

        VX1=VX2BOUND
        VY1=VY2BOUND
        VZ1=VZ2BOUND

        AX1=AX2
        AY1=AY2
        AZ1=AZ2

        BX1=BX2
        BY1=BY2
        BZ1=BZ2

        IF (ISNORDER.EQ.0) THEN

          X2B=X1+VX1*DT
          Y2B=Y1+VY1*DT
          Z2B=Z1+VZ1*DT

        ELSE

          CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX1,BY1,BZ1,DT,
     &      X2B,Y2B,Z2B,
     &      VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM,GAMMA,ICHARGE,BMOVECUT,
     &      IUSTEP,IENELOSS,DGAMMA)

        ENDIF

        CALL MYBFELD(X2B,Y2B,Z2B,BX2,BY2,BZ2,AX2,AY2,AZ2)

        CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DTIM,
     &    X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,
     &    IUSTEP,IENELOSS,DGAMMA)

        IF (IENELOSS.NE.0) THEN
          DGAMSUM=DGAMSUM+VXSIGN*DGAMMA
          IF (ABS(DGAMSUM).GT.GAMMA*1.0D-8) THEN
            GAMMA=GAMMA+DGAMSUM
            DGAMSUM=0.0D0
          ENDIF
          BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
          VN=SQRT(VX2*VX2+VY2*VY2+VZ2*VZ2)
          VX2=VX2/VN*CLIGHT1*BETA
          VY2=VY2/VN*CLIGHT1*BETA
          VZ2=VZ2/VN*CLIGHT1*BETA
        ENDIF

        DIST1=X1-X2SAV
        DIST2=X2-X2SAV

        DDT=DTIM*DABS(DIST1)/(DABS(DIST1)+DABS(DIST2))
        DDT2=DDT*BSHIFT

        IF(ISNORDER.EQ.0) THEN

          X2B=X1+VX1*DDT2
          Y2B=Y1+VY1*DDT2
          Z2B=Z1+VZ1*DDT2

        ELSE

          CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX1,BY1,BZ1,DDT2,
     &      X2B,Y2B,Z2B,
     &      VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM,GAMMA,ICHARGE,BMOVECUT,
     &      IUSTEP,IENELOSS,DGAMMA)

        ENDIF

        CALL MYBFELD(X2B,Y2B,Z2B,BX2,BY2,BZ2,AX2,AY2,AZ2)

        CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDT,
     &    X2INT,Y2INT,Z2INT,VX2INT,VY2INT,VZ2INT,
     &    VXPINT,VYPINT,VZPINT,GAMMA,ICHARGE,BMOVECUT,
     &    IUSTEP,IENELOSS,DGAMMA)

        DISTI=X2INT-X2SAV

        DDDT=DDT
        IF (DIST1.NE.0.) DDDT=DDT-DDT*DISTI/DABS(DIST1)

        CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDDT,
     &    X3INT,Y3INT,Z3INT,VX3INT,VY3INT,VZ3INT,
     &    VXPINT,VYPINT,VZPINT,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

        DISTI=X3INT-X2SAV
        DDDDT=DDDT
        IF (DIST1.NE.0.) DDDDT=DDDT-DDDT*DISTI/DABS(DIST1)
        DDDDT2=BSHIFT*DDDDT

        CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDDDT,
     &    X2,Y2,Z2,VX2,VY2,VZ2,
     &    VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,
     &    IUSTEP,IENELOSS,DGAMMA)

        IF (IENELOSS.NE.0) THEN
          DGAMSUM=DGAMSUM+VXSIGN*DGAMMA
          IF (ABS(DGAMSUM).GT.GAMMA*1.0D-8) THEN
            GAMMA=GAMMA+DGAMSUM
            DGAMSUM=0.0D0
          ENDIF
          BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
          VN=SQRT(VX2*VX2+VY2*VY2+VZ2*VZ2)
          VX2=VX2/VN*CLIGHT1*BETA
          VY2=VY2/VN*CLIGHT1*BETA
          VZ2=VZ2/VN*CLIGHT1*BETA
        ENDIF

        X1=X1SAV
        Y1=Y1SAV
        Z1=Z1SAV

        T1=T1SAV

        VX1=VX1SAV
        VY1=VY1SAV
        VZ1=VZ1SAV

        AX1=AX1SAV
        AY1=AY1SAV
        AZ1=AZ1SAV

        BX1=BX1SAV
        BY1=BY1SAV
        BZ1=BZ1SAV

        DTIM=DTIM0
        DT=DTIM*BSHIFT

        IROI1=IROI1+1
        IROI2=IROI1+1

      ELSE ! BOUNDARY CROSSED

        IF (IENELOSS.NE.0) THEN
          DGAMSUM=DGAMSUM+VXSIGN*DGAMMA
          IF (ABS(DGAMSUM).GT.GAMMA*1.0D-8) THEN
            GAMMA=GAMMA+DGAMSUM
            DGAMSUM=0.0D0
          ENDIF
          BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
          VN=SQRT(VX2*VX2+VY2*VY2+VZ2*VZ2)
          VX2=VX2/VN*CLIGHT1*BETA
          VY2=VY2/VN*CLIGHT1*BETA
          VZ2=VZ2/VN*CLIGHT1*BETA
        ENDIF

      ENDIF ! BOUNDARY CROSSED

C--- BOUNDARY CROSSING }

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

        CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX1,BY1,BZ1,DDT2,
     &    X2B,Y2B,Z2B,
     &    VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM,GAMMA,ICHARGE,BMOVECUT,
     &    IUSTEP,IENELOSS,DGAMMA)

      ENDIF

      CALL MYBFELD(X2B,Y2B,Z2B,BX2,BY2,BZ2,AX2,AY2,AZ2)

      CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDT,
     &  X2INT,Y2INT,Z2INT,VX2INT,VY2INT,VZ2INT,
     &  VXPINT,VYPINT,VZPINT,GAMMA,ICHARGE,BMOVECUT,
     &  IUSTEP,IENELOSS,DGAMMA)

      DISTI=(X2INT-XF0)*EWSFX+(Y2INT-YF0)*EWSFY+(Z2INT-ZF0)*EWSFZ
      DDDT=DDT
      IF (DIST1.NE.0.) DDDT=DDT-DDT*DISTI/DABS(DIST1)

      CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDDT,
     &  X3INT,Y3INT,Z3INT,VX3INT,VY3INT,VZ3INT,
     &  VXPINT,VYPINT,VZPINT,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

      DISTI=(X3INT-XF0)*EWSFX+(Y3INT-YF0)*EWSFY+(Z3INT-ZF0)*EWSFZ
      DDDDT=DDDT
      IF (DIST1.NE.0.) DDDDT=DDDT-DDDT*DISTI/DABS(DIST1)
      DDDDT2=BSHIFT*DDDDT

      CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DDDDT,
     &  X2,Y2,Z2,VX2,VY2,VZ2,
     &  VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,
     &  IUSTEP,IENELOSS,DGAMMA)

      IF (IENELOSS.NE.0) THEN
        DGAMSUM=DGAMSUM+VXSIGN*DGAMMA
        IF (ABS(DGAMSUM).GT.GAMMA*1.0D-8) THEN
          GAMMA=GAMMA+DGAMSUM
          DGAMSUM=0.0D0
        ENDIF
        BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
        VN=SQRT(VX2*VX2+VY2*VY2+VZ2*VZ2)
        VX2=VX2/VN*CLIGHT1*BETA
        VY2=VY2/VN*CLIGHT1*BETA
        VZ2=VZ2/VN*CLIGHT1*BETA
      ENDIF

      T2=T1+DDDDT

      X1=XOLD
      Y1=YOLD
      Z1=ZOLD
      VX1=VXOLD
      VY1=VYOLD
      VZ1=VZOLD

      IF (IENELOSS.NE.0) THEN
        DGAMSUM=DGAMSUM+VXSIGN*DGAMMA
        GAMMA=GAMMA+DGAMSUM
        DGAMSUM=0.0D0
        BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
        VN=SQRT(VX2BOUND*VX2BOUND+VY2BOUND*VY2BOUND+VZ2BOUND*VZ2BOUND)
        VX2BOUND=VX2BOUND/VN*CLIGHT1*BETA
        VY2BOUND=VY2BOUND/VN*CLIGHT1*BETA
        VZ2BOUND=VZ2BOUND/VN*CLIGHT1*BETA
      ENDIF

      GAMMAL=GAMMA0-GAMMA

      RETURN
      END
