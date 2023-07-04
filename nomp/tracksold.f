*CMZ :  3.05/06 17/07/2018  11.15.17  by  Michael Scheer
*CMZ :  3.05/03 16/05/2018  16.01.42  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  11.18.26  by  Michael Scheer
*CMZ :  2.61/02 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  10.48.09  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  11.27.13  by  Michael Scheer
*CMZ :  2.16/08 20/10/2000  11.44.50  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.12/02 15/06/99  10.22.14  by  Michael Scheer
*CMZ :  2.12/01 14/06/99  15.14.15  by  Michael Scheer
*CMZ :  2.10/01 18/02/99  11.14.28  by  Michael Scheer
*CMZ :  1.04/03 11/12/98  11.42.46  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  09.52.03  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.43.20  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.32  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE TRACKSOLD(ISOUR)
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
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEND.

         INTEGER ISOUR,IZAEHL,IC,JC

          DOUBLE PRECISION X1,Y1,Z1,VX1,VY1,VZ1,BX1,BY1,BZ1,VXP,VYP,VZP
          DOUBLE PRECISION X2,Y2,Z2,VX2,VY2,VZ2,BX2,BY2,BZ2,BSQ,BS
          DOUBLE PRECISION X2B,Y2B,Z2B,AX2,AY2,AZ2
          DOUBLE PRECISION VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM
          DOUBLE PRECISION XENDSOU,DTIM,DT2,ECDUM

C--- START OF TRACKING

         X1=SOURCEA(1,1,ISOUR)
         Y1=SOURCEA(2,1,ISOUR)
         Z1=SOURCEA(3,1,ISOUR)

         VX1=SOURCEA(1,2,ISOUR)
         VY1=SOURCEA(2,2,ISOUR)
         VZ1=SOURCEA(3,2,ISOUR)

CERR10.1292    BX1=SOURCEA(1,3,ISOUR)
CERR10.1292    BY1=SOURCEA(2,3,ISOUR)
CERR10.1292    BZ1=SOURCEA(3,3,ISOUR)

          BX1=SOURCEA(1,4,ISOUR)
          BY1=SOURCEA(2,4,ISOUR)
          BZ1=SOURCEA(3,4,ISOUR)

CV2------------------------------------------------------
          IF (ISOUR.NE.ISOURO) THEN
CV2------------------------------------------------------
         ECSOUR(1,ISOUR)=0.0
         ECSOUR(2,ISOUR)=0.0
         ECSOUR(3,ISOUR)=0.0
         ECSOUR(4,ISOUR)=0.0
         ECMAX(   ISOUR)=-1.D30

CV2------------------------------------------------------
          IZTOT(ISOUR)=0
          ENDIF !ISOUR
CV2------------------------------------------------------
C--- STEP SIZE ACCORDING TO SOURCE LENGTH AND NUMBER OF STEPS

         XENDSOU=SOURCEE(1,1,ISOUR)    !FINAL X

C- ATTENTION: STEP SIZE IS DIFFERENT FOR EACH SOURCE POINT !!

         DTIM=(XENDSOU-X1)/NLPOI/CLIGHT1
         DT2=DTIM/2.D0

C- CHECK NUMBER OF STEPS

         IF (NLPOI/(XENDSOU-X1).LT.MYINUM) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'*** WARNING SR TRACKSOLD ***'
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'STEP SIZE FOR SOURCE POINT IS LOWER THAN STEP'
         WRITE(LUNGFO,*)'SIZE FOR TRAJECTORY!'
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'INCREASE NLPOI OR BE AWARE OF STRANGE RESULTS!'
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)'*** WARNING SR TRACKSOLD ***'
         WRITE(6,*)
         WRITE(6,*)'STEP SIZE FOR SOURCE POINT IS LOWER THAN STEP'
         WRITE(6,*)'SIZE FOR TRAJECTORY!'
         WRITE(6,*)
         WRITE(6,*)'INCREASE NLPOI OR BE AWARE OF STRANGE RESULTS!'
         WRITE(6,*)
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

1000     IZAEHL=IZAEHL+1

         IF (IZAEHL.GT.NDWSOU) THEN

         WRITE(LUNGFO,*)'*** ERROR IN SR TRACKSOLD ***'
         WRITE(LUNGFO,*)
     &'TOO MANY STEPS, INCREASE PARAMETER NBADDP IN SOURCE.CMN'
         WRITE(6,*)'*** ERROR IN SR TRACKSOLD ***'
         WRITE(6,*)
     &'TOO MANY STEPS, INCREASE PARAMETER NBADDP IN SOURCE.CMN'
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

         CALL BMOVE(X1,Y1,Z1,VX1,VY1,VZ1,BX1,BY1,BZ1,DT2,
     &     X2B,Y2B,Z2B,
     &     VXDUM,VYDUM,VZDUM,VXPDUM,VYPDUM,VZPDUM,DMYGAMMA,ICHARGE,BMOVECUT,IUSTEP)

      ENDIF

         CALL MYBFELD(X2B,Y2B,Z2B,BX2,BY2,BZ2,AX2,AY2,AZ2)

         CALL BMOVE(X1,Y1,Z1,VX1,VY1,VZ1,BX2,BY2,BZ2,DTIM,
     &     X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,DMYGAMMA,ICHARGE,BMOVECUT,IUSTEP)

C- STORE POINT

          WSOU(1,1,IZAEHL)=X2
          WSOU(2,1,IZAEHL)=Y2
          WSOU(3,1,IZAEHL)=Z2
        wsou(1,4,IZAEHL)=DTIM

          WSOU(1,2,IZAEHL)=VX2
          WSOU(2,2,IZAEHL)=VY2
          WSOU(3,2,IZAEHL)=VZ2

          WSOU(1,3,IZAEHL)=VXP
          WSOU(2,3,IZAEHL)=VYP
          WSOU(3,3,IZAEHL)=VZP

          BSQ=BX2**2+BY2**2+BZ2**2
          BS=DSQRT(BSQ)
          ECSOUR(1,ISOUR)=ECSOUR(1,ISOUR)+BS
          ECSOUR(4,ISOUR)=ECSOUR(4,ISOUR)+DSIGN(BS,BY2)
          ECSOUR(3,ISOUR)=ECSOUR(3,ISOUR)+BSQ
          IF (ECMAX(ISOUR).LT.DSQRT(BSQ)) ECMAX(ISOUR)=BS

C--- END OF LOOP

      IF (X2.LT.XENDSOU)  GOTO 1000

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
