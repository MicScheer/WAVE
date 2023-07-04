*CMZ :  3.07/00 05/03/2019  13.36.15  by  Michael Scheer
*CMZ :  3.02/00 19/09/2014  10.44.41  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.67/02 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.67/01 15/03/2012  17.06.56  by  Michael Scheer
*CMZ :  2.66/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.65/03 02/10/2009  13.09.16  by  Michael Scheer
*CMZ :  2.65/02 28/09/2009  12.44.41  by  Michael Scheer
*CMZ :  2.64/01 20/08/2009  11.37.22  by  Michael Scheer
*CMZ :  2.63/05 12/08/2009  08.49.28  by  Michael Scheer
*CMZ :  2.63/03 02/05/2008  14.41.00  by  Michael Scheer
*CMZ :  2.61/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  10.48.09  by  Michael Scheer
*CMZ :  2.52/02 08/07/2004  10.00.17  by  Michael Scheer
*CMZ :  2.34/07 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.31/01 24/04/2001  17.59.02  by  Michael Scheer
*CMZ :  2.30/02 12/04/2001  19.11.52  by  Michael Scheer
*CMZ :  2.30/01 12/04/2001  14.52.49  by  Michael Scheer
*CMZ :  2.20/12 11/04/2001  16.39.04  by  Michael Scheer
*CMZ :  2.20/11 11/04/2001  11.00.41  by  Michael Scheer
*CMZ :  2.20/10 10/04/2001  12.05.35  by  Michael Scheer
*CMZ :  2.20/09 03/04/2001  10.08.21  by  Michael Scheer
*CMZ :  2.15/01 30/03/2001  17.09.31  by  Michael Scheer
*CMZ :  2.20/08 18/03/2001  21.50.06  by  Michael Scheer
*CMZ :  2.20/07 18/03/2001  17.08.58  by  Michael Scheer
*CMZ :  2.20/06 15/03/2001  17.20.08  by  Michael Scheer
*CMZ :  2.20/05 15/03/2001  16.57.30  by  Michael Scheer
*CMZ :  2.20/04 09/03/2001  16.47.40  by  Michael Scheer
*CMZ :  2.20/03 23/02/2001  15.04.13  by  Michael Scheer
*CMZ :  2.20/02 21/02/2001  11.30.46  by  Michael Scheer
*CMZ :  2.20/01 20/02/2001  14.18.37  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE TRASOU(ISOUR)
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

C--- EVALUATE INTEGRALES FOR A SINGLE SOURCE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,primkin.
      include 'primkin.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEND.

      DOUBLE PRECISION X0,X1,X2,Y1,Y2,Z1,Z2,XENDSOU,TLEN
     &  ,T,DT,DT2,DT0,DT01,DTIM01,VXP,VYP,VZP
      DOUBLE PRECISION VX1,VY1,VZ1
      DOUBLE PRECISION VX2,VY2,VZ2,BX2,BY2,BZ2
      DOUBLE PRECISION C1,GAMMA0,GAMMA,DGAMSUM,VN,BETA,DGAMMA

      DOUBLE PRECISION X2B,Y2B,Z2B,AX2D,AY2D,AZ2D

      INTEGER ISOUR,IZAEHL,ICAL,LSTEP,MCOO

      DATA ICAL/0/
      DATA MCOO/0/

      IF (ICAL.EQ.0) THEN
        DTIM01=1.D0/DTIM0
        C1=1.D0/CLIGHT1
      ENDIF !ICAL

      GAMMA0=DMYGAMMA
      GAMMA=GAMMA0
      DGAMSUM=0.0D0

      X1=SOURCEAO(1,1,ISOUR)
      Y1=SOURCEAO(2,1,ISOUR)
      Z1=SOURCEAO(3,1,ISOUR)

      VX1=SOURCEAO(1,2,ISOUR)
      VY1=SOURCEAO(2,2,ISOUR)
      VZ1=SOURCEAO(3,2,ISOUR)

      XENDSOU=SOURCEEO(1,1,ISOUR)    !FINAL X
      X0=X1

      TLEN=SOURCET(2,ISOUR)-SOURCET(1,ISOUR)
c 20140918      DT0=TLEN/FLOAT(NCO-1)/FLOAT(NSOURCE)
      DT0=TLEN/FLOAT(nlpoi-1)/float(nsource)
      DTMCO=DT0
      DT01=1.D0/DT0

      MCO=TLEN/DT0+1

      IF (MCOO.EQ.0) THEN

        ALLOCATE(DWT(MCO))
        ALLOCATE(DWX(MCO))
        ALLOCATE(DWX2P(MCO))
        ALLOCATE(DWB(MCO))
        ALLOCATE(DWB2P(MCO))
        ALLOCATE(DWY(MCO))
        ALLOCATE(DWY2P(MCO))
        ALLOCATE(DWZ(MCO))
        ALLOCATE(DWZ2P(MCO))
        ALLOCATE(TRAGAM(MCO))

      ELSE IF (MCOO.LT.MCO) THEN

        DEALLOCATE(DWT)
        DEALLOCATE(DWX)
        DEALLOCATE(DWX2P)
        DEALLOCATE(DWB)
        DEALLOCATE(DWB2P)
        DEALLOCATE(DWY)
        DEALLOCATE(DWY2P)
        DEALLOCATE(DWZ)
        DEALLOCATE(DWZ2P)
        DEALLOCATE(TRAGAM)

        ALLOCATE(DWT(MCO))
        ALLOCATE(DWX(MCO))
        ALLOCATE(DWX2P(MCO))
        ALLOCATE(DWB(MCO))
        ALLOCATE(DWB2P(MCO))
        ALLOCATE(DWY(MCO))
        ALLOCATE(DWY2P(MCO))
        ALLOCATE(DWZ(MCO))
        ALLOCATE(DWZ2P(MCO))
        ALLOCATE(TRAGAM(MCO))

      ENDIF

      T=0.D0

      DT=DT0
      DT2=DT/2.D0

      X2=X1
      Y2=Y1
      Z2=Z1

      VX2=VX1
      VY2=VY1
      VZ2=VZ1

C LOOP OVER TIME STEPS

      IZAEHL=1

      DWT(IZAEHL)=T
      DWX(IZAEHL)=X1
      DWY(IZAEHL)=Y1
      DWZ(IZAEHL)=Z1

C GET MAGNETIC FIELD {

      X2B=X1+VX1*DT2
      Y2B=Y1+VY1*DT2
      Z2B=Z1+VZ1*DT2
      CALL MYBFELD(X2B,Y2B,Z2B,BX2,BY2,BZ2,AX2D,AY2D,AZ2D)

C GET MAGNETIC FIELD }

      DWB(IZAEHL)=SIGN(SQRT(BX2*BX2+BY2*BY2+BZ2*BZ2),BY2)
      TRAGAM(IZAEHL)=GAMMA

      LSTEP=0
1000  IZAEHL=IZAEHL+1

      IF (IZAEHL.EQ.MCO-1) THEN
        DT=(XENDSOU-X2)/VX2
        LSTEP=1
      ENDIF

      X1=X2
      Y1=Y2
      Z1=Z2

      VX1=VX2
      VY1=VY2
      VZ1=VZ2

      T=T+DT

C GET MAGNETIC FIELD {

      X2B=X1+VX1*DT2
      Y2B=Y1+VY1*DT2
      Z2B=Z1+VZ1*DT2
      CALL MYBFELD(X2B,Y2B,Z2B,BX2,BY2,BZ2,AX2D,AY2D,AZ2D)

C GET MAGNETIC FIELD }

C MOVE ONE STEP {

      CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DT,
     &  X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

      IF (IENELOSS.NE.0) THEN
        DGAMSUM=DGAMSUM+DGAMMA
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

C MOVE ONE STEP }

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION }

C--- END OF LOOP OVER TIME STEPS

      IF (IENELOSS.NE.0) THEN
        DGAMSUM=DGAMSUM+DGAMMA
        GAMMA=GAMMA+DGAMSUM
        DGAMSUM=0.0D0
        BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))
        VN=SQRT(VX2*VX2+VY2*VY2+VZ2*VZ2)
        VX2=VX2/VN*CLIGHT1*BETA
        VY2=VY2/VN*CLIGHT1*BETA
        VZ2=VZ2/VN*CLIGHT1*BETA
      ENDIF

      DWT(IZAEHL)=T
      DWX(IZAEHL)=X2
      DWB(IZAEHL)=SIGN(SQRT(BX2*BX2+BY2*BY2+BZ2*BZ2),BY2)
      DWY(IZAEHL)=Y2
      DWZ(IZAEHL)=Z2
      TRAGAM(IZAEHL)=GAMMA

      IF (IZAEHL.LT.MCO.AND.LSTEP.EQ.0)  GOTO 1000
      MCO=IZAEHL

      CALL util_spline_coef_F90(DWT,DWX,MCO,-9999.0d0,-9999.0d0,DWX2P)
      CALL util_spline_coef_F90(DWT,DWB,MCO,-9999.0d0,-9999.0d0,DWB2P)
      CALL util_spline_coef_F90(DWT,DWY,MCO,-9999.0d0,-9999.0d0,DWY2P)
      CALL util_spline_coef_F90(DWT,DWZ,MCO,-9999.0d0,-9999.0d0,DWZ2P)

      MCOO=MCO

      RETURN
      END
