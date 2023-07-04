*CMZ :  4.00/13 07/11/2021  15.08.37  by  Michael Scheer
*CMZ :  4.00/04 29/07/2019  16.45.12  by  Michael Scheer
*CMZ :  3.05/28 19/12/2018  10.06.30  by  Michael Scheer
*CMZ :  3.05/03 16/05/2018  16.01.42  by  Michael Scheer
*CMZ :  3.05/01 07/05/2018  16.35.20  by  Michael Scheer
*CMZ :  3.03/02 07/12/2015  17.18.10  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.66/18 30/11/2010  14.43.20  by  Michael Scheer
*CMZ :  2.65/02 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.54/04 14/09/2009  15.19.42  by  Michael Scheer
*CMZ :  2.54/03 18/04/2005  08.54.26  by  Michael Scheer
*CMZ :  2.50/03 10/05/2004  14.40.44  by  Michael Scheer
*CMZ :  2.50/00 28/04/2004  12.39.19  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  15.58.48  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  2.14/02 26/04/2000  16.42.15  by  Michael Scheer
*CMZ :  2.12/02 15/06/99  10.22.13  by  Michael Scheer
*CMZ :  2.11/01 21/05/99  10.11.03  by  Michael Scheer
*CMZ :  2.10/01 04/03/99  10.05.14  by  Michael Scheer
*CMZ :  2.02/00 05/02/99  14.51.22  by  Michael Scheer
*CMZ :  2.01/00 19/01/99  11.16.08  by  Michael Scheer
*CMZ :  2.00/00 04/01/99  14.07.05  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.43.04  by  Michael Scheer
*CMZ :  1.03/01 26/01/98  17.10.38  by  Michael Scheer
*CMZ :  1.03/00 16/01/98  13.56.25  by  Michael Scheer
*CMZ :  1.02/03 13/01/98  17.27.39  by  Michael Scheer
*CMZ : 00.02/05 04/03/97  13.43.34  by  Michael Scheer
*CMZ : 00.02/01 12/12/96  16.16.58  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.12.46  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.53.31  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.45  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE REARG(ISOUR,IOBSV)
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

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,reargf90u.
      include 'reargf90u.cmn'
*KEND.


C--- CALCULATES REAL PARTS OF INTEGRALS FOR A GIVEN OBSERVATION POINT
C    THE TIME STARTS BY ZERO, I.E. EACH SOURCE HAS ITS OWN TIME BASE
C    SOURCES ARE ADDED INCOHERENTLY

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,colli.
      include 'colli.cmn'
*KEND.

      INTEGER ISOUR,IOBSV,IPOI,IPOIR,INSIDE,INOLD,ICHANGE
      INTEGER IWARN,IRANGE,IAPERT,IWARNBET1N

      DOUBLE PRECISION C1,DT,BET1N,BET1NO
      DOUBLE PRECISION X,Y,Z,T,BX,BY,BZ,BPX,BPY,BPZ,ZP,YP,OPANG,TOPANG
      DOUBLE PRECISION XOB,YOB,ZOB,RX,RY,RZ
      DOUBLE PRECISION RNX,RNY,RNZ,RNBX,RNBY,RNBZ,R,R1,PX,PY,PZ,DOM1,DOM2
      DOUBLE PRECISION APX,APY,APZ,APERW2,APERH2,RARGV2(3)
      double precision br2,rnr2,br4,rnr4,b3

      DOUBLE PRECISION R0
      DATA R0/1.D30/

      DATA IWARN/0/
      DATA IWARNBET1N/1/ ! 7.12.2015, since warning is boring

      integer, save :: ksouro=-1

      save iwarn,iwarnbet1n

      IF (ISOUR.NE.ksouro) IWARN=0

      C1=1.D0/CLIGHT1

C--- LOOP OVER TIME INTERVALS

CV2--------------------------------------------------------------
CORR  IF (ISOUR.NE.ksouro) T=0.0d0
C      R0=OBSV(1,IOBSV)-SOURCEEO(1,1,ISOUR)   !C150793
      R0=OBSV(1,1)-SOURCEAO(1,1,ISOUR)   !C210599
      T=TBUFF(IOBSV)
CV2--------------------------------------------------------------

CV2   T=0.0d0

      OPANG=DABS(WGWINFC/DMYGAMMA)
      TOPANG=DTAN(OPANG)
      IPOIR=0
      INSIDE=0
      INOLD=0
      ICHANGE=0

      XOB=OBSV(1,IOBSV)
      YOB=OBSV(2,IOBSV)
      ZOB=OBSV(3,IOBSV)

      DO IPOI=1,IPOISOU(ISOUR)

        X=WSOU(1,1,IPOI)
        Y=WSOU(2,1,IPOI)
        Z=WSOU(3,1,IPOI)

        DT=WSOU(1,4,IPOI)

C DT FOR STEP FROM X(I-1) TO X(I) IS WSOU(1,4,I)
C HENCE TIME OF X(I) SHOULD BE:

        T=T+DT

        IF (T/DT.GT.1.D7) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN REARG ***'
          WRITE(LUNGFO,*)
     &      'T/DT.GT.1.D7 (TIME INTERVALL OF TRACKING STEP TOO SMALL)'
          WRITE(LUNGFO,*)'CHECK PROBLEM USING DEBUGGER'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN REARG ***'
          WRITE(6,*)
     &      'T/DT.GT.1.D7 (TIME INTERVALL OF TRACKING STEP TOO SMALL)'
          WRITE(6,*)'CHECK PROBLEM USING DEBUGGER'
          WRITE(6,*)
          WRITE(6,*)
        ENDIF    !T/DT

        BX=WSOU(1,2,IPOI)*C1
        BY=WSOU(2,2,IPOI)*C1
        BZ=WSOU(3,2,IPOI)*C1

        BPX=WSOU(1,3,IPOI)*C1
        BPY=WSOU(2,3,IPOI)*C1
        BPZ=WSOU(3,3,IPOI)*C1

        RX=XOB-X
        RY=YOB-Y
        RZ=ZOB-Z

        R=DSQRT(RX*RX+RY*RY+RZ*RZ)
        R1=1.D0/R

        RNX=RX*R1
        RNY=RY*R1
        RNZ=RZ*R1

C--- THE DISTANCE R IS INTRODUCED HERE EXPLICITLY (S. PROGRAM OF CHAOEN WANG

        BET1N=(1.0D0-BX*RNX)-BY*RNY-BZ*RNZ

c 20090928{
        br2=by**2+bz**2
        rnr2=rny**2+rnz**2
        b3=dmybeta**3
        br4=br2**2
        rnr4=rnr2**2

        if(br2.lt.1.0d-4.and.rnr2.lt.1.0d-4) then
          bet1n=
     &      1.0d0/(1+dmybeta)/dmygamma**2
     &      +dmybeta*(rnr2/2.0d0
     &      +rnr4/8.0d0)
     &      +(br2/2.0d0
     &      -br2*rnr2/4.0d0
     &      -br2*rnr4/16.0d0)/dmybeta
     &      +b3*br4*(1.0d0/8.0d0
     &      -rnr2/16.0d0
     &      -rnr4/64.0d0)
     &      -by*rny
     &      -bz*rnz
        endif
c }20090928
        DOM1=1.D0/(R*BET1N*BET1N)

        RNBX=RNX-BX
        RNBY=RNY-BY
        RNBZ=RNZ-BZ

        PX=(RNBY*BPZ-RNBZ*BPY)
        PY=(RNBZ*BPX-RNBX*BPZ)
        PZ=(RNBX*BPY-RNBY*BPX)

        ZP=BZ/BX
        YP=BY/BX

        IF(
     &      Z+(ZP+TOPANG)/(1.D0-ZP*TOPANG)*RX.GE.ZOB
     &      .AND.
     &      Z+(ZP-TOPANG)/(1.D0+ZP*TOPANG)*RX.LE.ZOB
     &      .AND.
     &      Y+(YP+TOPANG)/(1.D0-YP*TOPANG)*RX.GE.YOB
     &      .AND.
     &      Y+(YP-TOPANG)/(1.D0+YP*TOPANG)*RX.LE.YOB)
     &      THEN

          IPOIR=IPOIR+1

          IF(IPOIR.GT.NDARGU) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN REARG ***'
            WRITE(LUNGFO,*)'INCREASE PARAMETER NDARGUP IN CMPARA.CMN'
            WRITE(LUNGFO,*)
            WRITE(6,*) '*** ERROR IN REARG ***'
            STOP
          ENDIF

          INSIDE=1

          APX=(APERX-X)
          APZ=Z+RZ/RX*APX
          APY=Y+RY/RX*APX
          APERW2=APERWID/2.D0
          APERH2=APERHIG/2.D0

          IF (X.GE.XIANF.AND.X.LE.XIEND) THEN
            IRANGE=1
          ELSE
            IRANGE=0
          ENDIF

          IF (
     &        APZ.LE.(APERZ+APERW2)
     &        .AND.
     &        APZ.GE.(APERZ-APERW2)
     &        .AND.
     &        APY.LE.(APERY+APERH2)
     &        .AND.
     &        APY.GE.(APERY-APERH2)
     &        ) THEN
            IAPERT=1
          ELSE
            IAPERT=0
          ENDIF


          IF (IRANGE.EQ.1.AND.IAPERT.EQ.1) THEN


            IF (IVELOFIELD.EQ.0) THEN
              DOM2=CLIGHT1*DOM1*R1/(DMYGAMMA*DMYGAMMA)
              REARGUM(1,IPOIR)=(RNY*PZ-RNZ*PY)*DOM1+(RNX-BX)*DOM2
              REARGUM(2,IPOIR)=(RNZ*PX-RNX*PZ)*DOM1+(RNY-BY)*DOM2
              REARGUM(3,IPOIR)=(RNX*PY-RNY*PX)*DOM1+(RNZ-BZ)*DOM2
            ELSE IF (IVELOFIELD.EQ.1) THEN
              REARGUM(1,IPOIR)=(RNY*PZ-RNZ*PY)*DOM1
              REARGUM(2,IPOIR)=(RNZ*PX-RNX*PZ)*DOM1
              REARGUM(3,IPOIR)=(RNX*PY-RNY*PX)*DOM1
            ELSE IF (IVELOFIELD.EQ.2) THEN !approximation, see souintana
              REARGUM(1,IPOIR)=RNBX*R1
              REARGUM(2,IPOIR)=RNBY*R1
              REARGUM(3,IPOIR)=RNBZ*R1
            ELSE IF (IVELOFIELD.LT.0) THEN
              DOM2=CLIGHT1*DOM1*R1/(DMYGAMMA*DMYGAMMA)
              REARGUM(1,IPOIR)=(RNX-BX)*DOM2
              REARGUM(2,IPOIR)=(RNY-BY)*DOM2
              REARGUM(3,IPOIR)=(RNZ-BZ)*DOM2
            ELSE  !IVELOFIELD
              PRINT*, '*** ERROR IN REARG: BAD VALUE OF IVELOFIELD'
              STOP
            ENDIF !IVELOFIELD

C150793             REARGUM(4,IPOIR)=T+R*C1
            REARGUM(4,IPOIR)=T+(R-R0)*C1

            IF (IVELOFIELD.NE.2) THEN

              REARGUM(5,IPOIR)=(
     &          REARGUM(1,IPOIR)*REARGUM(1,IPOIR)
     &          +REARGUM(2,IPOIR)*REARGUM(2,IPOIR)
     &          +REARGUM(3,IPOIR)*REARGUM(3,IPOIR)
     &          )/BET1N

            ELSE

              DOM2=CLIGHT1*DOM1*R1/(DMYGAMMA*DMYGAMMA)
              RARGV2(1)=(RNY*PZ-RNZ*PY)*DOM1+(RNX-BX)*DOM2
              RARGV2(2)=(RNZ*PX-RNX*PZ)*DOM1+(RNY-BY)*DOM2
              RARGV2(3)=(RNX*PY-RNY*PX)*DOM1+(RNZ-BZ)*DOM2

              REARGUM(5,IPOIR)=(
     &          RARGV2(1)*RARGV2(1)
     &          +RARGV2(2)*RARGV2(2)
     &          +RARGV2(3)*RARGV2(3)
     &          )/BET1N

            ENDIF

            REARGUM(6,IPOIR)=DT
            REARGUM(7,IPOIR)=BET1N
            REARGUM(8,IPOIR)=-CLIGHT1*R1
     &        *((BX*BX-(BX*RNX+BY*RNY+BZ*RNZ)**2)+BY*BY+BZ*BZ)
     &        -(RNX*BPX+RNY*BPY+RNZ*BPZ)/2.0D0

          ELSE

            IF (IAPERT.EQ.0.AND.IWARN.EQ.0) THEN

              WRITE(LUNGFO,*)'*** WARNING REARG:'
              WRITE(LUNGFO,*)
     &          'APERTURE AFFECTS SPECTRUM CALCULATION FOR SOURCE'
     &          ,ISOUR
              WRITE(LUNGFO,*)
              WRITE(6,*)'*** WARNING REARG:'
              WRITE(6,*)
     &          'APERTURE AFFECTS SPECTRUM CALCULATION FOR SOURCE'
     &          ,ISOUR
              WRITE(6,*)
              IWARN=1

            ENDIF    !IAPERT

            REARGUM(1,IPOIR)=0.0d0
            REARGUM(2,IPOIR)=0.0d0
            REARGUM(3,IPOIR)=0.0d0

            REARGUM(4,IPOIR)=T+(R-R0)*C1
            REARGUM(5,IPOIR)=0.0d0
            REARGUM(6,IPOIR)=DT

            REARGUM(7,IPOIR)=0.0d0
            REARGUM(8,IPOIR)=0.0d0

          ENDIF    !IRANGE, IAPERT

          IF (IPOIR.NE.1) THEN
            IF ((R-R0)*C1/T.GT.1.D7) THEN

              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** WARNING REARG ***'
              WRITE(LUNGFO,*)
     &          '(R-R0)*C1/T.GT.1.D7'
              WRITE(LUNGFO,*)'CHECK PROBLEM USING DEBUGGER'
              WRITE(LUNGFO,*)
              WRITE(6,*)
              WRITE(6,*)
              WRITE(6,*)'*** WARNING REARG ***'
              WRITE(6,*)
     &          '(R-R0)*C1/T.GT.1.D7'
              WRITE(6,*)'CHECK PROBLEM USING DEBUGGER'
              WRITE(6,*)
              WRITE(6,*)

            ENDIF    !T/(R-R0)*C1
          ENDIF    !IPOIR.NE.1

        ELSE

          INSIDE=0

        ENDIF   !INSIDE WGWINFC CONE

        IF(INOLD.NE.INSIDE) THEN
          ICHANGE=ICHANGE+1
        ENDIF

        IF(ICHANGE.GT.2) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN REARG ***'
          WRITE(LUNGFO,*)'SOURCE POINT NOT COHERENT FOR OBSERVATION POINT'
          WRITE(LUNGFO,*)'X-START AND  X-END OF SOURCE:'
          WRITE(LUNGFO,*)SOURCEAO(1,1,ISOUR),SOURCEEO(1,1,ISOUR)
          WRITE(LUNGFO,*)'X,Y,Z OF OBSERVATION POINT:'
          WRITE(LUNGFO,*)XOB,YOB,ZOB
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'CHECK OBSERVATION POINTS, OR INCREASE MYINUM,'
          WRITE(LUNGFO,*)'OR CHANGE WGWINFC, OR CHANGE COLLIMATOR ...'
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN REARG ***'
          WRITE(6,*)'SOURCE POINT NOT COHERENT FOR OBSERVATION POINT'
          WRITE(6,*)'X-START AND  X-END OF SOURCE:'
          WRITE(6,*)SOURCEAO(1,1,ISOUR),SOURCEEO(1,1,ISOUR)
          WRITE(6,*)'X,Y,Z OF OBSERVATION POINT:'
          WRITE(6,*)XOB,YOB,ZOB
          WRITE(6,*)
          WRITE(6,*)'CHECK OBSERVATION POINTS OR INCREASE MYINUM,'
          WRITE(6,*)'OR CHANGE WGWINFC, OR CHANGE COLLIMATOR ...'
          WRITE(6,*)
          STOP
        ENDIF

        INOLD=INSIDE

      ENDDO !LOOP OVER TIME INTERVALLS

      IARGUM=IPOIR
CV2------------------------------------------------------------
      IF (ISOUR.NE.ksouro) NARGUM(IOBSV,ISOUR)=0
      NARGUM(IOBSV,ISOUR)=NARGUM(IOBSV,ISOUR)+IARGUM
      TBUFF(IOBSV)=T
CV2------------------------------------------------------------
CV2   NARGUM(IOBSV,ISOUR)=IARGUM

      ksouro=isour

      RETURN
      END
