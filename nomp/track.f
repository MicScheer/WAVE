*CMZ :  4.00/11 22/11/2020  13.40.46  by  Michael Scheer
*CMZ :  3.07/00 04/03/2019  18.43.21  by  Michael Scheer
*CMZ :  3.06/00 21/01/2019  12.37.35  by  Michael Scheer
*CMZ :  3.05/04 28/06/2018  09.36.14  by  Michael Scheer
*CMZ :  3.03/04 18/12/2017  11.38.05  by  Michael Scheer
*CMZ :  3.02/00 27/08/2014  16.23.17  by  Michael Scheer
*CMZ :  3.01/06 20/06/2014  16.35.06  by  Michael Scheer
*CMZ :  3.01/00 26/06/2013  16.44.13  by  Michael Scheer
*CMZ :  3.00/01 28/03/2013  09.49.39  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.69/00 25/10/2012  15.53.04  by  Michael Scheer
*CMZ :  2.67/05 15/05/2012  14.54.44  by  Michael Scheer
*CMZ :  2.67/02 03/05/2012  14.14.09  by  Michael Scheer
*CMZ :  2.66/13 03/06/2010  08.58.45  by  Michael Scheer
*CMZ :  2.66/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.65/02 28/09/2009  12.44.41  by  Michael Scheer
*CMZ :  2.64/00 14/08/2009  14.46.01  by  Michael Scheer
*CMZ :  2.63/05 12/08/2009  14.40.36  by  Michael Scheer
*CMZ :  2.61/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.59/01 24/01/2007  13.53.39  by  Michael Scheer
*CMZ :  2.57/05 10/01/2007  13.38.27  by  Michael Scheer
*CMZ :  2.57/02 30/11/2005  11.58.01  by  Michael Scheer
*CMZ :  2.54/07 16/06/2005  11.40.58  by  Michael Scheer
*CMZ :  2.53/04 09/02/2005  11.05.58  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  10.48.09  by  Michael Scheer
*CMZ :  2.52/08 14/10/2004  14.52.40  by  Michael Scheer
*CMZ :  2.50/00 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.47/12 03/07/2003  09.48.54  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  11.09.32  by  Michael Scheer
*CMZ :  2.16/08 31/10/2000  17.33.47  by  Michael Scheer
*CMZ :  2.16/06 27/08/2000  19.50.51  by  Michael Scheer
*CMZ :  2.16/05 28/07/2000  17.03.57  by  Michael Scheer
*CMZ :  2.16/04 19/07/2000  10.39.46  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.13/10 23/03/2000  13.09.47  by  Michael Scheer
*CMZ :  2.13/11 22/03/2000  15.55.53  by  Michael Scheer
*CMZ :  2.12/01 08/06/99  15.38.34  by  Michael Scheer
*CMZ :  2.12/00 02/06/99  13.56.19  by  Michael Scheer
*CMZ :  2.11/00 11/05/99  15.47.07  by  Michael Scheer
*CMZ :  2.10/01 17/02/99  14.03.33  by  Michael Scheer
*CMZ :  1.04/00 11/12/98  10.40.37  by  Michael Scheer
*CMZ :  1.00/00 11/07/97  15.11.58  by  Michael Scheer
*CMZ : 00.01/10 16/07/96  14.53.50  by  Michael Scheer
*CMZ : 00.01/09 11/04/96  17.41.50  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  09.51.12  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.41.55  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.31  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE TRACK(X1,Y1,Z1,VX1,VY1,VZ1,
     &  XF0,YF0,ZF0,EWSFX,EWSFY,EWSFZ,
     &  X2,Y2,Z2,VX2,VY2,VZ2,DTIM,BSHIFT,GAMMA0,GAMMAL)
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

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEND.

C---     INPUT:
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

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,optic.
      include 'optic.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,primkin.
      include 'primkin.cmn'
*KEEP,photon.
      include 'photon.cmn'
*KEEP,ustep.
      include 'ustep.cmn'
*KEND.

      INTEGER ICAL,IZAEHL,I,J,K,I1,I2,IROI1,IROI2,istat,nfit,ifit1,ifit2

      DOUBLE PRECISION X1,Y1,Z1,VX1,VY1,VZ1,X2,Y2,Z2,VX2,VY2,VZ2
     &  ,DTIM,BSHIFT,X2B,Y2B,Z2B,BX1,BY1,BZ1,BX2,BY2,BZ2
     &  ,GAMMA0,GAMMA,DT,VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM
     &  ,X2INT,Y2INT,Z2INT,DDT,DDDT,DDT2
     &  ,VX2INT,VY2INT,VZ2INT,VXPINT,VYPINT,VZPINT
     &  ,VXP,VYP,VZP,XOLD,YOLD,ZOLD,VXOLD,VYOLD,VZOLD
     &  ,X3INT,Y3INT,Z3INT,VX3INT,VY3INT,VZ3INT,DDDDT,DDDDT2
     &  ,EWSFX,EWSFY,EWSFZ,XF0,YF0,ZF0,DIST1,DIST2,DISTI
     &  ,AX1,AY1,AZ1,AX2,AY2,AZ2,DXYZ,PDUM,BBMAX,phdum
     &  ,X2BOUND,Y2BOUND,Z2BOUND,VX2BOUND,VY2BOUND,VZ2BOUND
     &  ,X1SAV,Y1SAV,Z1SAV,VX1SAV,VY1SAV,VZ1SAV,BX1SAV,BY1SAV,BZ1SAV
     &  ,AX1SAV,AY1SAV,AZ1SAV,T1SAV,X2SAV,BPER(3),VN
     &  ,DGAMMA,GAMMAL,BETA,DGAMSUM,dumf

      DOUBLE PRECISION T1,T2,DW
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xfit,yfit,efit
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: WVP

      DATA ICAL/0/

      nutrack=-1
      nustep=-1

      GAMMA=GAMMA0
      DGAMMA=0.0D0
      DGAMSUM=0.0D0

      VXP=0.0D0
      VYP=0.0D0
      VZP=0.0D0

      IF (ICAL.EQ.0) THEN

        XMX=-1.D30
        YMX=-1.D30
        ZMX=-1.D30
        BXMX=-1.D30
        BYMX=-1.D30
        BZMX=-1.D30
        XMXAE=-1.D30
        YMXAE=-1.D30
        ZMXAE=-1.D30
        BXMXAE=-1.D30
        BYMXAE=-1.D30
        BZMXAE=-1.D30
        PHIMX=-1.D30
        BBMAX=-1.D30

        XMN=1.D30
        YMN=1.D30
        ZMN=1.D30
        XMNAE=1.D30
        YMNAE=1.D30
        ZMNAE=1.D30
        BXMN=1.D30
        BYMN=1.D30
        BZMN=1.D30
        BXMNAE=1.D30
        BYMNAE=1.D30
        BZMNAE=1.D30
        PHIMN=1.D30

      ENDIF   !ICAL

      IF (ICAL.EQ.0) THEN
        ALLOCATE(WSXYZ(3,NSTEPMX))
        ALLOCATE(WVXYZ(3,NSTEPMX))
        ALLOCATE(WBXYZ(3,NSTEPMX))
        ALLOCATE(WAXYZ(3,NSTEPMX))
        ALLOCATE(WGAMMA(3,0:NSTEPMX))
        WGAMMA(2,1)=0.0D0
        ALLOCATE(WVP(3,NSTEPMX))
        ALLOCATE(WTIM0(NSTEPMX))
        ALLOCATE(HTRA2(NSTEPMX))
        ALLOCATE(WTRA2(NSTEPMX))
      ENDIF !ICAL

c        GAMMA22=1.0D0/(GAMMA*GAMMA*2.D0)

      IF (ICAL.EQ.0) THEN

        BINT0=0.0
        BINT1X=0.0
        BINT1Y=0.0
        BINT1Z=0.0
        BINT2X=0.0
        BINT2Y=0.0
        BINT2Z=0.0
        BINT3Y=0.0
        BINT3YA=0.0
        PINT=0.0
        B2INTY=0.0
        WTRA2I=0.0D0

      ENDIF !ICAL

      WTRA2IC=0.0D0

      DTIM0=DTIM

c        PDUM=1.D9/2./PI1*CGAM1*CLIGHT1**2/1.D18*DMYCUR*DMYENERGY**2
      PDUM=0.5d-9/PI1*CGAM1*DMYCUR*(CLIGHT1*EMASSG1)**2

      DT=DTIM*BSHIFT

      T1=0.0

      CALL MYBFELD(X1,Y1,Z1,BX1,BY1,BZ1,AX1,AY1,AZ1)

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

      VN=SQRT(VX1*VX1+VY1*VY1+VZ1*VZ1)

      VX2=VX1
      VY2=VY1
      VZ2=VZ1

      BX2=BX1
      BY2=BY1
      BZ2=BZ1

      IZAEHL=0

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

1000  IZAEHL=IZAEHL+1

      nustep=izaehl

      IF (IZAEHL.GT.NWMAX) THEN

        WRITE(LUNGFO,*)
     &    '*** WARNING SR TRACK: TOO MANY STEPS, INCREASE PARAMETER NWMAXP IN CMPARA.CMN ***'
        WRITE(6,*)
     &    '*** WARNING SR TRACK: TOO MANY STEPS, INCREASE PARAMETER NWMAXP IN CMPARA.CMN ***'

        IZAEHL=IZAEHL-1
        XSTOP=WSXYZ(1,IZAEHL)

        WRITE(LUNGFO,*)'     NUMBER OF STEPS DONE:',IZAEHL
        WRITE(6,*)'     NUMBER OF STEPS DONE:',IZAEHL
        WRITE(LUNGFO,*)
     &    '*** XSTOP overwritten:',XSTOP,' ***'
        WRITE(6,*)
     &    '*** XSTOP overwritten:',XSTOP,' ***'

        GOTO 1001

      ENDIF !(IZAEHL.GT.NWMAX)

      IF(ICAL.EQ.0) THEN

C-- INTEGRALS OF B-FIELD

        DXYZ=DSQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1)+(Z2-Z1)*(Z2-Z1)) !18.12.91
        BINT0=BINT0+DXYZ

        IF (X2.NE.X1) THEN
          DW=0.5D0*(X2-X1)*(((Y2-Y1)/(X2-X1))**2+((Z2-Z1)/(X2-X1))**2)
        ELSE IF (Y2.NE.Y1.OR.Z2.NE.Z1) THEN
          DW=1.D30
        ENDIF

        WTRA2I=WTRA2I+DW

        IF (IZAEHL.EQ.1) THEN
          WTRA2(IZAEHL)=DW
          HTRA2(IZAEHL)=DW+(X2-X1)*(1.0D0/(GAMMA*GAMMA*2.D0))
        ELSE
          WTRA2(IZAEHL)=WTRA2(IZAEHL-1)+DW
          HTRA2(IZAEHL)=HTRA2(IZAEHL-1)+DW+(X2-X1)*(1.0D0/(GAMMA*GAMMA*2.D0))
        ENDIF

        BINT1X=BINT1X+BX2*DXYZ
        BINT1Y=BINT1Y+BY2*DXYZ
        BINT1Z=BINT1Z+BZ2*DXYZ
        BINT2X=BINT2X+BINT1X*DXYZ
        BINT2Y=BINT2Y+BINT1Y*DXYZ
        BINT2Z=BINT2Z+BINT1Z*DXYZ
        BINT3YA=BINT3YA+DABS(BY2*BY2*BY2)*DXYZ
        BINT3Y=BINT3Y+BY2*BY2*BY2*DXYZ

        BPER(1)=(VY2*BZ2-VZ2*BY2)/VN
        BPER(2)=(VZ2*BX2-VX2*BZ2)/VN
        BPER(3)=(VX2*BY2-VY2*BX2)/VN

        PINT=PINT+
     &    PDUM*GAMMA**2*
     &    (BPER(1)**2+BPER(2)**2+BPER(3)**2)*DXYZ !LEISTUNG/BAHN

C-- PEAK-VALUES

        PHIX=DATAN(VZ2/VX2)

        IF (PHIMX.LT.PHIX) PHIMX=PHIX
        IF (PHIMN.GT.PHIX) PHIMN=PHIX

        IF (XMX.LT.X2) XMX=X2
        IF (XMN.GT.X2) XMN=X2
        IF (YMX.LT.Y2) YMX=Y2
        IF (YMN.GT.Y2) YMN=Y2
        IF (ZMX.LT.Z2) ZMX=Z2
        IF (ZMN.GT.Z2) ZMN=Z2

        IF (BXMX.LT.BX2) BXMX=BX2
        IF (BXMN.GT.BX2) BXMN=BX2
        IF (BYMX.LT.BY2) BYMX=BY2
        IF (BYMN.GT.BY2) BYMN=BY2
        IF (BZMX.LT.BZ2) BZMX=BZ2
        IF (BZMN.GT.BZ2) BZMN=BZ2

        IF (XIANF.LE.X2 .AND. XIEND.GE.X2) THEN
          IF (XMXAE.LT.X2) XMXAE=X2
          IF (XMNAE.GT.X2) XMNAE=X2
          IF (YMXAE.LT.Y2) YMXAE=Y2
          IF (YMNAE.GT.Y2) YMNAE=Y2
          IF (ZMXAE.LT.Z2) ZMXAE=Z2
          IF (ZMNAE.GT.Z2) ZMNAE=Z2

          IF (BXMXAE.LT.BX2) BXMXAE=BX2
          IF (BXMNAE.GT.BX2) BXMNAE=BX2
          IF (BYMXAE.LT.BY2) BYMXAE=BY2
          IF (BYMNAE.GT.BY2) BYMNAE=BY2
          IF (BZMXAE.LT.BZ2) BZMXAE=BZ2
          IF (BZMNAE.GT.BZ2) BZMNAE=BZ2
        ENDIF

        BRMS=BX2*BX2+BY2*BY2+BZ2*BZ2

        IF (BRMS.GT.BBMAX) BBMAX=BRMS

        WSXYZ(1,IZAEHL)=X2
        WSXYZ(2,IZAEHL)=Y2
        WSXYZ(3,IZAEHL)=Z2

        WVXYZ(1,IZAEHL)=VX2
        WVXYZ(2,IZAEHL)=VY2
        WVXYZ(3,IZAEHL)=VZ2

        WBXYZ(1,IZAEHL)=BX2
        WBXYZ(2,IZAEHL)=BY2
        WBXYZ(3,IZAEHL)=BZ2

        WAXYZ(1,IZAEHL)=AX2
        WAXYZ(2,IZAEHL)=AY2
        WAXYZ(3,IZAEHL)=AZ2

        WVP(1,IZAEHL)=VXP
        WVP(2,IZAEHL)=VYP
        WVP(3,IZAEHL)=VZP

        WTIM0(IZAEHL)=T2

        WGAMMA(1,IZAEHL)=GAMMA !absolute value
        WGAMMA(2,IZAEHL)=WGAMMA(2,IZAEHL-1)+DGAMMA !accumulated loss
        WGAMMA(3,IZAEHL)=DGAMMA !loss per step

C        WRITE (LUNTR,*)X2,Y2,Z2,VX2,VY2,VZ2,BX2,BY2,BZ2,AX2,AY2,AZ2

      ENDIF   !ICAL

C        WRITE (LUNTR,*)X2,Y2,Z2,VX2,VY2,VZ2,BX2,BY2,BZ2,AX2,AY2,AZ2

      X1=X2
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

      ITRACK=IZAEHL

      XTRACK=X1
      YTRACK=Y1
      ZTRACK=Z1

      VXTRACK=VX1
      VYTRACK=VY1
      VZTRACK=VZ1

      IF (ISNORDER.EQ.0) THEN

        X2B=X1+VX1*DT
        Y2B=Y1+VY1*DT
        Z2B=Z1+VZ1*DT

      ELSE

        CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX1,BY1,BZ1,DT,
     &    X2B,Y2B,Z2B,
     &    VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM,GAMMA,ICHARGE,BMOVECUT,
     &    IUSTEP,IENELOSS,DGAMMA)

      ENDIF

      CALL MYBFELD(X2B,Y2B,Z2B,BX2,BY2,BZ2,AX2,AY2,AZ2)

      CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DTIM,
     &  X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,
     &  IUSTEP,IENELOSS,DGAMMA)

      DXYZ=DSQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1)+(Z2-Z1)*(Z2-Z1))

      IF (X2.NE.X1) THEN
        DW=0.5D0*(X2-X1)*(((Y2-Y1)/(X2-X1))**2+((Z2-Z1)/(X2-X1))**2)
      ELSE IF (Y2.NE.Y1.OR.Z2.NE.Z1) THEN
        DW=1.0D30
      ENDIF

      WTRA2IC=WTRA2IC+DW

C--- BOUNDARY CROSSING {

      IF (X2.GT.ROIX(IROI2)) THEN

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
     &    VXPINT,VYPINT,VZPINT,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

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

        IF (IENELOSS.gt.0) THEN

          DGAMSUM=DGAMSUM+DGAMMA

          IF (ABS(DGAMSUM).GT.GAMMA*1.0D-8) THEN

            phdum=emasse1*DSQRT((gamma-1.0d0)*(gamma+1.0d0))/vn

            GAMMA=GAMMA+DGAMSUM
            DGAMSUM=0.0D0
            BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
            VN=SQRT(VX2bound*VX2bound+VY2bound*VY2bound+VZ2bound*VZ2bound)

            VX2bound=VX2bound/VN*CLIGHT1*BETA
            VY2bound=VY2bound/VN*CLIGHT1*BETA
            VZ2bound=VZ2bound/VN*CLIGHT1*BETA

          ENDIF ! summe

        ELSE IF (IENELOSS.lt.0) THEN

          phdum=emasse1*DSQRT((gamma-1.0d0)*(gamma+1.0d0))/vn

          GAMMA=GAMMA+DGAMMA
          BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
          VN=SQRT(VX2bound*VX2bound+VY2bound*VY2bound+VZ2bound*VZ2bound)

cerror ?? 21.1.2019
c          vx2bound=(phdum*vx2bound-pphoton(1))/(emasse1*gamma)*clight1
c          vy2bound=(phdum*vy2bound-pphoton(1))/(emasse1*gamma)*clight1
c          vz2bound=(phdum*vz2bound-pphoton(1))/(emasse1*gamma)*clight1

cerror ?? 21.1.2019
          vx2bound=(phdum*vx2bound-dpphoton(1))/(emasse1*gamma)*clight1
          vy2bound=(phdum*vy2bound-dpphoton(2))/(emasse1*gamma)*clight1
          vz2bound=(phdum*vz2bound-dpphoton(3))/(emasse1*gamma)*clight1

        ENDIF !ieneloss.ne.0

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

          DGAMSUM=DGAMSUM+DGAMMA

          IF (ABS(DGAMSUM).GT.GAMMA*1.0D-8) THEN

            phdum=emasse1*DSQRT((gamma-1.0d0)*(gamma+1.0d0))/vn

            GAMMA=GAMMA+DGAMSUM
            DGAMSUM=0.0D0
            BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
            VN=SQRT(VX2*VX2+VY2*VY2+VZ2*VZ2)

            IF (IENELOSS.eq.-1) THEN
              vx2=(phdum*vx2-pphoton(1))/(emasse1*gamma)*clight1
              vy2=(phdum*vy2-pphoton(1))/(emasse1*gamma)*clight1
              vz2=(phdum*vz2-pphoton(1))/(emasse1*gamma)*clight1
            else
              VX2=VX2/VN*CLIGHT1*BETA
              VY2=VY2/VN*CLIGHT1*BETA
              VZ2=VZ2/VN*CLIGHT1*BETA
            endif

          ENDIF ! summe

        ENDIF !ieneloss.ne.0

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
     &    VXPINT,VYPINT,VZPINT,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

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
     &    VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

        IF (IENELOSS.NE.0) THEN

          DGAMSUM=DGAMSUM+DGAMMA

          IF (ABS(DGAMSUM).GT.GAMMA*1.0D-8) THEN

            phdum=emasse1*DSQRT((gamma-1.0d0)*(gamma+1.0d0))/vn

            GAMMA=GAMMA+DGAMSUM
            DGAMSUM=0.0D0
            BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
            VN=SQRT(VX2*VX2+VY2*VY2+VZ2*VZ2)

            IF (IENELOSS.eq.-1) THEN
              vx2=(phdum*vx2-pphoton(1))/(emasse1*gamma)*clight1
              vy2=(phdum*vy2-pphoton(1))/(emasse1*gamma)*clight1
              vz2=(phdum*vz2-pphoton(1))/(emasse1*gamma)*clight1
            else
              VX2=VX2/VN*CLIGHT1*BETA
              VY2=VY2/VN*CLIGHT1*BETA
              VZ2=VZ2/VN*CLIGHT1*BETA
            endif

          ENDIF ! summe

        ENDIF !ieneloss.ne.0

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

          DGAMSUM=DGAMSUM+DGAMMA

          IF (ABS(DGAMSUM).GT.GAMMA*1.0D-8) THEN

            phdum=emasse1*DSQRT((gamma-1.0d0)*(gamma+1.0d0))/vn

            GAMMA=GAMMA+DGAMSUM
            DGAMSUM=0.0D0
            BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
            VN=SQRT(VX2*VX2+VY2*VY2+VZ2*VZ2)

            IF (IENELOSS.eq.-1) THEN
              vx2=(phdum*vx2-pphoton(1))/(emasse1*gamma)*clight1
              vy2=(phdum*vy2-pphoton(1))/(emasse1*gamma)*clight1
              vz2=(phdum*vz2-pphoton(1))/(emasse1*gamma)*clight1
            else
              VX2=VX2/VN*CLIGHT1*BETA
              VY2=VY2/VN*CLIGHT1*BETA
              VZ2=VZ2/VN*CLIGHT1*BETA
            endif

          ENDIF ! summe

        ENDIF !ieneloss.ne.0

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

        CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX1,BY1,BZ1,DDT2,X2B,Y2B,Z2B,
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
     &  VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

      IF (IENELOSS.NE.0) THEN

        DGAMSUM=DGAMSUM+DGAMMA

        IF (ABS(DGAMSUM).GT.GAMMA*1.0D-8) THEN

          phdum=emasse1*DSQRT((gamma-1.0d0)*(gamma+1.0d0))/vn

          GAMMA=GAMMA+DGAMSUM
          DGAMSUM=0.0D0
          BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
          VN=SQRT(VX2*VX2+VY2*VY2+VZ2*VZ2)

          IF (IENELOSS.eq.-1) THEN
            vx2=(phdum*vx2-pphoton(1))/(emasse1*gamma)*clight1
            vy2=(phdum*vy2-pphoton(1))/(emasse1*gamma)*clight1
            vz2=(phdum*vz2-pphoton(1))/(emasse1*gamma)*clight1
          else
            VX2=VX2/VN*CLIGHT1*BETA
            VY2=VY2/VN*CLIGHT1*BETA
            VZ2=VZ2/VN*CLIGHT1*BETA
          endif

        ENDIF ! summe

      ENDIF !ieneloss.ne.0

      T2=T1+DDDDT

C-- INTEGRALS OF B-FIELD

      DXYZ=DSQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1)+(Z2-Z1)*(Z2-Z1)) !18.12.91

      IF (ICAL.EQ.0) THEN

        BINT0=BINT0+DXYZ

        IF (X2.NE.X1) THEN
          DW=0.5D0*(X2-X1)*(((Y2-Y1)/(X2-X1))**2+((Z2-Z1)/(X2-X1))**2)
        ELSE IF (Y2.NE.Y1.OR.Z2.NE.Z1) THEN
          DW=1.D30
        ENDIF

        WTRA2I=WTRA2I+DW

        IF (IZAEHL.EQ.1) THEN
          HTRA2(IZAEHL)=DW+(X2-X1)*(1.0D0/(GAMMA*GAMMA*2.D0))
          WTRA2(IZAEHL)=DW
        ELSE
          WTRA2(IZAEHL)=WTRA2(IZAEHL-1)+DW
          HTRA2(IZAEHL)=HTRA2(IZAEHL-1)+DW+(X2-X1)*(1.0D0/(GAMMA*GAMMA*2.D0))
        ENDIF

        BINT1X=BINT1X+BX2*DXYZ
        BINT1Y=BINT1Y+BY2*DXYZ
        BINT1Z=BINT1Z+BZ2*DXYZ
        BINT2X=BINT2X+BINT1X*DXYZ
        BINT2Y=BINT2Y+BINT1Y*DXYZ
        BINT2Z=BINT2Z+BINT1Z*DXYZ
        BINT3YA=BINT3YA+DABS(BY2*BY2*BY2)*DXYZ
        BINT3Y=BINT3Y+BY2*BY2*BY2*DXYZ
cerror 26.6.2013        PINT=PINT+PDUM*(BX2*BX2+BY2*BY2+BZ2*BZ2)*DXYZ !LEISTUNG/BAHN
        PINT=PINT+PDUM*GAMMA**2*(BX2*BX2+BY2*BY2+BZ2*BZ2)*DXYZ !LEISTUNG/BAHN
        B2INTY=PINT/PDUM/GAMMA**2

C-- PEAK-VALUES

        PHIX=DATAN(VZ2/VX2)

        IF (PHIMX.LT.PHIX) PHIMX=PHIX
        IF (PHIMN.GT.PHIX) PHIMN=PHIX

        IF (XMX.LT.X2) XMX=X2
        IF (XMN.GT.X2) XMN=X2
        IF (YMX.LT.Y2) YMX=Y2
        IF (YMN.GT.Y2) YMN=Y2
        IF (ZMX.LT.Z2) ZMX=Z2
        IF (ZMN.GT.Z2) ZMN=Z2

        IF (BXMX.LT.BX2) BXMX=BX2
        IF (BXMN.GT.BX2) BXMN=BX2
        IF (BYMX.LT.BY2) BYMX=BY2
        IF (BYMN.GT.BY2) BYMN=BY2
        IF (BZMX.LT.BZ2) BZMX=BZ2
        IF (BZMN.GT.BZ2) BZMN=BZ2

      ENDIF !ICAL

      IF (X2.NE.X1) THEN
        DW=0.5D0*(X2-X1)*(((Y2-Y1)/(X2-X1))**2+((Z2-Z1)/(X2-X1))**2)
      ELSE IF (Y2.NE.Y1.OR.Z2.NE.Z1) THEN
        DW=1.D30
      ENDIF

      WTRA2IC=WTRA2IC+DW

1001  X1=XOLD
      Y1=YOLD
      Z1=ZOLD
      VX1=VXOLD
      VY1=VYOLD
      VZ1=VZOLD

      IF (ICAL.EQ.0) THEN

        NCO=IZAEHL

        ALLOCATE(WTRA(3,5,NCO))

        DO I=1,IZAEHL-1

          WTRA(1,1,I)=WSXYZ(1,I)
          WTRA(2,1,I)=WSXYZ(2,I)
          WTRA(3,1,I)=WSXYZ(3,I)

          WTRA(1,2,I)=WVXYZ(1,I)
          WTRA(2,2,I)=WVXYZ(2,I)
          WTRA(3,2,I)=WVXYZ(3,I)

          WTRA(1,3,I)=WBXYZ(1,I)
          WTRA(2,3,I)=WBXYZ(2,I)
          WTRA(3,3,I)=WBXYZ(3,I)

          WTRA(1,4,I)=WAXYZ(1,I)
          WTRA(2,4,I)=WAXYZ(2,I)
          WTRA(3,4,I)=WAXYZ(3,I)

          WTRA(1,5,I)=WGAMMA(1,I)
          WTRA(2,5,I)=WGAMMA(2,I)
          WTRA(3,5,I)=WGAMMA(3,I)

        ENDDO !IZAEHL

        WTIM0(IZAEHL)=T2

        WTRA(1,1,IZAEHL)=X2
        WTRA(2,1,IZAEHL)=Y2
        WTRA(3,1,IZAEHL)=Z2

        WTRA(1,2,IZAEHL)=VX2
        WTRA(2,2,IZAEHL)=VY2
        WTRA(3,2,IZAEHL)=VZ2

        WTRA(1,3,IZAEHL)=BX2
        WTRA(2,3,IZAEHL)=BY2
        WTRA(3,3,IZAEHL)=BZ2

        WTRA(1,4,IZAEHL)=AX2
        WTRA(2,4,IZAEHL)=AY2
        WTRA(3,4,IZAEHL)=AZ2

        WTRA(1,5,IZAEHL)=GAMMA
        WTRA(2,5,IZAEHL)=WTRA(2,5,IZAEHL-1)+DGAMMA
        WTRA(3,5,IZAEHL)=DGAMMA

        WVP(1,IZAEHL)=VXP
        WVP(2,IZAEHL)=VYP
        WVP(3,IZAEHL)=VZP



        DEALLOCATE(WSXYZ)
        DEALLOCATE(WVXYZ)
        DEALLOCATE(WBXYZ)
        DEALLOCATE(WAXYZ)
        DEALLOCATE(WGAMMA)

        ALLOCATE(WSXYZ(3,NCO))
        ALLOCATE(WVXYZ(3,NCO))
        ALLOCATE(WBXYZ(3,NCO))
        ALLOCATE(WAXYZ(3,NCO))
        ALLOCATE(WVPXYZ(3,NCO))

        DO I=1,IZAEHL

          WSXYZ(1,I)=WTRA(1,1,I)
          WSXYZ(2,I)=WTRA(2,1,I)
          WSXYZ(3,I)=WTRA(3,1,I)

          WVXYZ(1,I)=WTRA(1,2,I)
          WVXYZ(2,I)=WTRA(2,2,I)
          WVXYZ(3,I)=WTRA(3,2,I)

          WBXYZ(1,I)=WTRA(1,3,I)
          WBXYZ(2,I)=WTRA(2,3,I)
          WBXYZ(3,I)=WTRA(3,3,I)

          WAXYZ(1,I)=WTRA(1,4,I)
          WAXYZ(2,I)=WTRA(2,4,I)
          WAXYZ(3,I)=WTRA(3,4,I)

          WVPXYZ(1,I)=WVP(1,I)
          WVPXYZ(2,I)=WVP(2,I)
          WVPXYZ(3,I)=WVP(3,I)

        ENDDO !IZAEHL

        DO I=1,IZAEHL
          WTRA(1,1,I)=WTIM0(I)
          WTRA(2,1,I)=HTRA2(I)
          WTRA(3,1,I)=WTRA2(I)
        ENDDO  !IZAEHL

        DEALLOCATE(WVP)
        DEALLOCATE(WTIM0)
        DEALLOCATE(HTRA2)
        DEALLOCATE(WTRA2)

        ALLOCATE(WTIM0(NCO))
        ALLOCATE(HTRA2(NCO))
        ALLOCATE(WTRA2(NCO))

        DO I=1,IZAEHL
          WTIM0(I)=WTRA(1,1,I)
          HTRA2(I)=WTRA(2,1,I)
          WTRA2(I)=WTRA(3,1,I)
        ENDDO  !IZAEHL

        DO I=1,IZAEHL

          WTRA(1,1,I)=WSXYZ(1,I)
          WTRA(2,1,I)=WSXYZ(2,I)
          WTRA(3,1,I)=WSXYZ(3,I)

          WTRA(1,2,I)=WVXYZ(1,I)
          WTRA(2,2,I)=WVXYZ(2,I)
          WTRA(3,2,I)=WVXYZ(3,I)

          WTRA(1,3,I)=WBXYZ(1,I)
          WTRA(2,3,I)=WBXYZ(2,I)
          WTRA(3,3,I)=WBXYZ(3,I)

          WTRA(1,4,I)=WAXYZ(1,I)
          WTRA(2,4,I)=WAXYZ(2,I)
          WTRA(3,4,I)=WAXYZ(3,I)

        ENDDO !IZAEHL
        DO I=1,nco
          if (wsxyz(1,i).ge.xianf.and.ifit1.eq.0) ifit1=i
          if (wsxyz(1,i).le.xiend) ifit2=i
        enddo
        nfit=ifit2-ifit1+1

        ALLOCATE(xfit(nfit))
        ALLOCATE(yfit(nfit))
        ALLOCATE(efit(nfit))

        DO I=ifit1,ifit2
          xfit(i-ifit1+1)=wsxyz(1,I)
          yfit(i-ifit1+1)=wsxyz(2,I)
          efit(i-ifit1+1)=0.0d0
        enddo

        yslopetr=0.0d0
        yoffstr=0.0d0
        call util_straight_line_fit(nfit,xfit,yfit,efit,yslopetr,yoffstr,
     &    dumf,yslopeetr,yoffsetr,istat)

        if (istat.ne.0) then
          write(lungfo,*)'*** Warning in TRACK: Fit of straight line failed ***'
          write(6,*)'*** Warning in TRACK: Fit of straight line failed ***'
          yslopetr=0.0d0
          yoffstr=0.0d0
        endif

        DO I=ifit1,ifit2
          yfit(i-ifit1+1)=wsxyz(3,I)
          efit(i-ifit1+1)=0.0d0
        enddo

        call util_straight_line_fit(nfit,xfit,yfit,efit,zslopetr,zoffstr,
     &    dumf,zslopeetr,zoffsetr,istat)
        if (istat.ne.0) then
          write(lungfo,*)'*** Warning in TRACK: Fit of straight line failed ***'
          write(6,*)'*** Warning in TRACK: Fit of straight line failed ***'
          zslopetr=0.0d0
          zoffstr=0.0d0
        endif

        DO I=ifit1,ifit2
          yfit(i-ifit1+1)=wvxyz(2,I)/wvxyz(1,i)
          efit(i-ifit1+1)=0.0d0
        enddo

        call util_straight_line_fit(nfit,xfit,yfit,efit,ypslopetr,ypoffstr,
     &    dumf,ypslopeetr,ypoffsetr,istat)
        if (istat.ne.0) then
          write(lungfo,*)'*** Warning in TRACK: Fit of straight line failed ***'
          write(6,*)'*** Warning in TRACK: Fit of straight line failed ***'
          ypslopetr=0.0d0
          yoffstr=0.0d0
        endif

        DO I=ifit1,ifit2
          yfit(i-ifit1+1)=wvxyz(3,I)/wvxyz(1,i)
          efit(i-ifit1+1)=0.0d0
        enddo

        call util_straight_line_fit(nfit,xfit,yfit,efit,zpslopetr,zpoffstr,
     &    dumf,zpslopeetr,zpoffsetr,istat)
        if (istat.ne.0) then
          write(lungfo,*)'*** Warning in TRACK: Fit of straight line failed ***'
          write(6,*)'*** Warning in TRACK: Fit of straight line failed ***'
          zpslopetr=0.0d0
          zpoffstr=0.0d0
        endif

        deALLOCATE(xfit)
        deALLOCATE(yfit)
        deALLOCATE(efit)


        I1=1
        I2=IZAEHL
        BRMS=0.D0
        ANGRMS=0.D0

        IF (BBMAX.GT.0.0) THEN
          DO I=1,IZAEHL
            BY1=
     &        +WBXYZ(1,I)*WBXYZ(1,I)
     &        +WBXYZ(2,I)*WBXYZ(2,I)
     &        +WBXYZ(3,I)*WBXYZ(3,I)
            IF (BY1.GE.0.01*BBMAX.AND.I1.EQ.1) I1=I
            IF (BY1.GE.0.01*BBMAX) I2=I
          ENDDO   !IZAEHL
        ENDIF  !(BBMAX.GT.0.0)

        DO I=I1,I2
          BY1=
     &      +WBXYZ(1,I)*WBXYZ(1,I)
     &      +WBXYZ(2,I)*WBXYZ(2,I)
     &      +WBXYZ(3,I)*WBXYZ(3,I)
          BRMS=BRMS+BY1
          ANGRMS=ANGRMS
     &      +(WVXYZ(3,I)*WVXYZ(3,I)+WVXYZ(2,I)*WVXYZ(2,I))
     &      /(WVXYZ(1,I)*WVXYZ(1,I))
        ENDDO !I1,I2

        BRMS=SQRT(BRMS/(I2-I1+1))
        ANGRMS=SQRT(ANGRMS/(I2-I1+1))

        IF(IWFILT0.NE.0) THEN

          OPEN(UNIT=LUNTR,FILE=FILETR,STATUS='NEW')

          IF(IWFILT0.EQ.1) THEN
            WRITE(LUNTR,*)'* ',ICODE
            WRITE(LUNTR,*)'* ',CODE
          ELSE IF(IWFILT0.GT.1) THEN
            WRITE(LUNTR,*)ICODE
            WRITE(LUNTR,*)CODE
          ENDIF

          DO K=1,IZAEHL,ABS(IWFILT0)
            WRITE(LUNTR,*)((WTRA(I,J,K),I=1,3),J=1,4)
          ENDDO

          CLOSE (LUNTR)

          WRITE(6,*) 'REFERENCE ORBIT WRITTEN TO'
          WRITE(6,*) 'FILE:',FILETR

        ENDIF !WFILT0

        ICAL=1

      ENDIF ! ICAL

      HTRA2I=WTRA2I+(XSTOP-XSTART)*(1.0D0/(GAMMA*GAMMA*2.D0))
      TTRA2I=HTRA2I/CLIGHT1

      IF (BXMX.EQ.0.D0.AND.BXMN.EQ.0.D0.AND.BZMX.EQ.0.D0.AND.BZMN.EQ.0.D0) THEN
        IBYONLY=1
      ELSE
        IBYONLY=0
      ENDIF

      if (ieneloss.ne.0) then

        phdum=emasse1*DSQRT((gamma-1.0d0)*(gamma+1.0d0))/vn

        GAMMA=GAMMA+DGAMSUM
        DGAMSUM=0.0D0
        BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
        VN=SQRT(VX2*VX2+VY2*VY2+VZ2*VZ2)

        IF (IENELOSS.eq.-1) THEN
          vx2=(phdum*vx2-pphoton(1))/(emasse1*gamma)*clight1
          vy2=(phdum*vy2-pphoton(1))/(emasse1*gamma)*clight1
          vz2=(phdum*vz2-pphoton(1))/(emasse1*gamma)*clight1
        else
          VX2=VX2/VN*CLIGHT1*BETA
          VY2=VY2/VN*CLIGHT1*BETA
          VZ2=VZ2/VN*CLIGHT1*BETA
        endif

      endif

      GAMMAL=GAMMA0-GAMMA

      tint=t2

      RETURN
      END
