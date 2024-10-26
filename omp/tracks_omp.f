*CMZ :  4.00/13 06/12/2021  12.45.17  by  Michael Scheer
*CMZ :  3.07/00 06/12/2021  12.43.51  by  Michael Scheer
*CMZ :  3.05/28 18/12/2018  13.41.17  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.12.29  by  Michael Scheer
*CMZ :  3.05/03 17/05/2018  14.49.31  by  Michael Scheer
*CMZ :  3.05/02 09/05/2018  13.42.42  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.70/04 21/12/2012  11.45.01  by  Michael Scheer
*CMZ :  2.69/02 05/11/2012  13.19.45  by  Michael Scheer
*CMZ :  2.69/01 31/10/2012  10.05.04  by  Michael Scheer
*CMZ :  2.69/00 30/10/2012  15.53.37  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  11.18.26  by  Michael Scheer
*CMZ :  2.66/18 01/12/2010  16.35.59  by  Michael Scheer
*CMZ :  2.66/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.65/02 28/09/2009  12.44.41  by  Michael Scheer
*CMZ :  2.64/01 18/08/2009  11.35.40  by  Michael Scheer
*CMZ :  2.63/05 17/08/2009  14.16.11  by  Michael Scheer
*CMZ :  2.61/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  10.59.41  by  Michael Scheer
*CMZ :  2.52/16 21/01/2005  17.09.20  by  Michael Scheer
*CMZ :  2.50/00 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.41/13 14/04/2004  13.21.28  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  11.26.17  by  Michael Scheer
*CMZ :  2.16/08 25/10/2000  12.21.53  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  16.16.11  by  Michael Scheer
*CMZ :  2.12/02 15/06/99  10.22.14  by  Michael Scheer
*CMZ :  2.12/00 02/06/99  13.58.22  by  Michael Scheer
*CMZ :  2.11/01 12/05/99  17.27.41  by  Michael Scheer
*CMZ :  2.11/00 12/05/99  12.12.31  by  Michael Scheer
*CMZ :  2.10/01 18/02/99  11.14.28  by  Michael Scheer
*CMZ :  1.04/03 11/12/98  11.42.46  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  09.52.03  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.43.20  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.32  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE tracks_omp(ISOUR)
*KEEP,GPLHINT.
*KEND.

      use bunchmod

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEND.
      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,optic.
      include 'optic.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,reargf90.
      include 'reargf90.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,primkin.
      include 'primkin.cmn'
*KEND.

      INTEGER ISOUR,IZAEHL,IC,JC,IWARN,iroi,kroi

      DOUBLE PRECISION X1,Y1,Z1,VX1,VY1,VZ1,BX1,BY1,BZ1,VXP,VYP,VZP
      DOUBLE PRECISION X2,Y2,Z2,VX2,VY2,VZ2,BX2,BY2,BZ2,BSQ,BS
      DOUBLE PRECISION X2B,Y2B,Z2B
      DOUBLE PRECISION VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM
      DOUBLE PRECISION XENDSOU,DTIM,DT2,ECDUM,dtim00,tlen,AX2D,AY2D,AZ2D
      DOUBLE PRECISION GAMMA,DGAMSUM,DGAMMA,BETA,VN
      integer lstep

      DATA IWARN/0/

C--- START OF TRACKING

      dgamsum=0.0d0 !? 6.12.2021

      if (ibunch.eq.0) then

        GAMMA=SOURCEG(1,1,ISOUR)
        BETA=DSQRT((1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA))

        X1=SOURCEA(1,1,ISOUR)
        Y1=SOURCEA(2,1,ISOUR)
        Z1=SOURCEA(3,1,ISOUR)

        VX1=SOURCEA(1,2,ISOUR)
        VY1=SOURCEA(2,2,ISOUR)
        VZ1=SOURCEA(3,2,ISOUR)

      else !ibunch

        gamma=egamma
        beta=dsqrt((1.0d0-1.0d0/gamma)*(1.d0+1.0d0/gamma))

        x1=xelec
        y1=yelec
        z1=zelec

        vx1=vxelec
        vy1=vyelec
        vz1=vzelec

      endif !ibunch

      BX1=SOURCEA(1,4,ISOUR)
      BY1=SOURCEA(2,4,ISOUR)
      BZ1=SOURCEA(3,4,ISOUR)

      lstep=0

      kroi=1
      do iroi=1,nroi
        if (x1.gt.roix(iroi)) kroi=iroi
      enddo

      IF (ISOUR.NE.ISOURO) THEN

        ECSOUR(1,ISOUR)=0.0
        ECSOUR(2,ISOUR)=0.0
        ECSOUR(3,ISOUR)=0.0
        ECSOUR(4,ISOUR)=0.0
        ECMAX(   ISOUR)=-1.D30

        IZTOT(ISOUR)=0
      ENDIF !ISOUR

C--- STEP SIZE ACCORDING TO SOURCE LENGTH AND NUMBER OF STEPS

      XENDSOU=SOURCEE(1,1,ISOUR)    !FINAL X

C- ATTENTION: STEP SIZE IS DIFFERENT FOR EACH SOURCE POINT !!

c      DTIM=(XENDSOU-X1)/(NLPOI-1)/CLIGHT1
      tlen=sourcet(2,isour)-sourcet(1,isour)
      dtim=tlen/(nlpoi-1) !in TRASOU wird hier durch NSOURCE dividiert!??
      dtim00=dtim
      DT2=DTIM/2.D0

C- CHECK NUMBER OF STEPS

      IF (IWARN.EQ.0) THEN
        IF (NLPOIO/(SOURCEEO(1,1,NSOURCE)-SOURCEAO(1,1,1)).LT.MYINUM) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING IN tracks_omp ***'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'STEP SIZE FOR SOURCE POINT IS LOWER THAN STEP'
          WRITE(LUNGFO,*)'SIZE FOR TRAJECTORY!'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'INCREASE NLPOI OR BE AWARE OF STRANGE RESULTS!'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING IN tracks_omp ***'
          WRITE(6,*)
          WRITE(6,*)'STEP SIZE FOR SOURCE POINT IS LOWER THAN STEP'
          WRITE(6,*)'SIZE FOR TRAJECTORY!'
          WRITE(6,*)
          WRITE(6,*)'INCREASE NLPOI OR BE AWARE OF STRANGE RESULTS!'
          WRITE(6,*)
        ENDIF
        IWARN=1
      ENDIF


      X2=X1
      Y2=Y1
      Z2=Z1

      VX2=VX1
      VY2=VY1
      VZ2=VZ1

      BX2=BX1
      BY2=BY1
      BZ2=BZ1

      IZAEHL=0 !LOOP COUNTER

C--- LOOP OVER STEPS

1000  IZAEHL=IZAEHL+1

      if (x2.gt.roix(kroi+1)) then
        kroi=kroi+1
        dtim=dtim00/roip(kroi)
      endif

      IF (LSTEP.EQ.1) THEN

        IF (X2.LE.XENDSOU) THEN

          DTIM=(XENDSOU-X2)/VX2/roip(kroi)
          DT2=DTIM/2.0D0

        ELSE

          DTIM=(XENDSOU-X1)/VX2/roip(kroi)
          DT2=DTIM/2.0D0

          X2=X1
          Y2=Y1
          Z2=Z1

          VX2=VX1
          VY2=VY1
          VZ2=VZ1

          BX2=BX1
          BY2=BY1
          BZ2=BZ1

          IZAEHL=IZAEHL-1

        ENDIF

      ENDIF !LSTEP.EQ.1

      IF (IZAEHL.GT.NDWSOU) THEN

        print*,IZAEHL,NDWSOU
        print*,X2

        WRITE(LUNGFO,*)'*** ERROR IN SR tracks_omp ***'
        WRITE(LUNGFO,*)
     &    'TOO MANY STEPS, INCREASE PARAMETER NBADDP IN SOURCE.CMN'
        WRITE(6,*)'*** ERROR IN SR tracks_omp ***'
        WRITE(6,*)
     &    'TOO MANY STEPS, INCREASE PARAMETER NBADDP IN SOURCE.CMN'
        STOP

      ENDIF

      X1=X2
      Y1=Y2
      Z1=Z2

      VX1=VX2
      VY1=VY2
      VZ1=VZ2

      IF (ISNORDER.EQ.0) THEN

        BX1=BX2
        BY1=BY2
        BZ1=BZ2

        X2B=X1+VX1*DT2
        Y2B=Y1+VY1*DT2
        Z2B=Z1+VZ1*DT2

      ELSE

        CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX1,BY1,BZ1,DT2,
     &    X2B,Y2B,Z2B,
     &    VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM,GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

      ENDIF

      CALL MYBFELD(X2B,Y2B,Z2B,BX2,BY2,BZ2,AX2D,AY2D,AZ2D)

C MYBFELD}

      CALL BMOVETAYL(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DTIM,
     &  X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,
     &  GAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

      IF (IENELOSS.NE.0.and.lstep.eq.0) THEN
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

C- STORE POINT

      WSOU(1,1,IZAEHL)=X2
      WSOU(2,1,IZAEHL)=Y2
      WSOU(3,1,IZAEHL)=Z2

      WSOU(1,2,IZAEHL)=VX2
      WSOU(2,2,IZAEHL)=VY2
      WSOU(3,2,IZAEHL)=VZ2

      WSOU(1,3,IZAEHL)=VXP
      WSOU(2,3,IZAEHL)=VYP
      WSOU(3,3,IZAEHL)=VZP

      wsou(1,4,IZAEHL)=DTIM
      wsou(2,4,IZAEHL)=BETA
      wsou(3,4,IZAEHL)=GAMMA

      wsou(1,5,IZAEHL)=bx2
      wsou(2,5,IZAEHL)=by2
      wsou(3,5,IZAEHL)=bz2

      BSQ=BX2*BX2+BY2*BY2+BZ2*BZ2
      BS=SQRT(BSQ)
      ECSOUR(1,ISOUR)=ECSOUR(1,ISOUR)+BS
      ECSOUR(4,ISOUR)=ECSOUR(4,ISOUR)+DSIGN(BS,BY2)
      ECSOUR(3,ISOUR)=ECSOUR(3,ISOUR)+BSQ

      IF (ECMAX(ISOUR).LT.DSQRT(BSQ)) ECMAX(ISOUR)=BS

C--- END OF LOOP

c      IF (X2.LT.XENDSOU)  GOTO 1000
      IF (X2.LT.XENDSOU-2.0d0*VX2*DTIM.AND.LSTEP.EQ.0)
     &  GOTO 1000

      IF (LSTEP.EQ.0) THEN
        LSTEP=1
        GOTO 1000
      ENDIF

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

C- STORE POINT

      WSOU(1,1,IZAEHL)=X2
      WSOU(2,1,IZAEHL)=Y2
      WSOU(3,1,IZAEHL)=Z2

      WSOU(1,2,IZAEHL)=VX2
      WSOU(2,2,IZAEHL)=VY2
      WSOU(3,2,IZAEHL)=VZ2

      WSOU(1,3,IZAEHL)=VXP
      WSOU(2,3,IZAEHL)=VYP
      WSOU(3,3,IZAEHL)=VZP

      wsou(1,4,IZAEHL)=DTIM
      wsou(2,4,IZAEHL)=BETA
      wsou(3,4,IZAEHL)=GAMMA

      wsou(1,5,IZAEHL)=bx2
      wsou(2,5,IZAEHL)=by2
      wsou(3,5,IZAEHL)=bz2

C- STORE NUMBER OF POINTS FOR INTEGRATION

      IPOISOU(ISOUR)=IZAEHL
CV2------------------------------------------------------
      IZTOT(ISOUR)=IZTOT(ISOUR)+IZAEHL
CERR101292      DO JC=1,3
      DO JC=1,2
        DO IC=1,3
          SOURCEA(IC,JC,ISOUR)=WSOU(IC,JC,IZAEHL)
        ENDDO
      ENDDO
      SOURCEA(1,4,ISOUR)=BX2
      SOURCEA(2,4,ISOUR)=BY2
      SOURCEA(3,4,ISOUR)=BZ2
CV2------------------------------------------------------


CV2------------------------------------------------------
      IF (NSADD.NE.0) THEN
CV2------------------------------------------------------

        ECSOUR(1,ISOUR)=ECSOUR(1,ISOUR)/IZTOT(ISOUR)
        ECSOUR(4,ISOUR)=ECSOUR(4,ISOUR)/IZTOT(ISOUR)
        ECSOUR(3,ISOUR)=ECSOUR(3,ISOUR)/IZTOT(ISOUR)
        ECDUM=ECSOUR(3,ISOUR)-ECSOUR(1,ISOUR)**2
        IF (ECDUM.LT.0.0) ECDUM=0.
        ECSOUR(3,ISOUR)=DSQRT(ECDUM)/ECSOUR(1,ISOUR)
        ECSOUR(2,ISOUR)=ECSOUR(1,ISOUR)*ecdipev1*DMYENERGY**2   !CRITICAL ENERGY
C260194  IF (IUNIT.NE.0)  ECSOUR(2,ISOUR)=ECSOUR(2,ISOUR)/WTOE1
        IF (IUNIT.NE.0)  ECSOUR(2,ISOUR)=WTOE1/ECSOUR(2,ISOUR)

CV2------------------------------------------------------
      ENDIF !NSADD
CV2------------------------------------------------------


      RETURN
      END
