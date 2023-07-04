*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.54/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.52/11 29/11/2004  16.35.20  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  1.03/06 06/08/98  18.01.07  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  14.48.44  by  Michael Scheer
*CMZ : 00.00/07 10/05/94  15.54.47  by  Michael Scheer
*CMZ : 00.00/05 29/04/94  19.47.02  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.07  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.02  by  Michael Scheer
*-- Author : Michael Scheer
C*****************************************************************
      SUBROUTINE BFELD( BX, BY, BZ,XX,Y,Z)
C*****************************************************************
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

C INPUT ORT X,Y,Z
C OUT B(X,Y,Z)

      IMPLICIT NONE

      EXTERNAL DCOSD,DSIND
      DOUBLE PRECISION DCOSD,DSIND

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,bfeld.
      include 'bfeld.cmn'
*KEND.

*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DOUBLE PRECISION BX,BY,BZ,XX,X,Y,Z,PELEV
      DOUBLE PRECISION XLBBY1,XLBBY2,XLBBY3,XLBBY4,XLBBY5,XLBBY6,XLBBY7
      DOUBLE PRECISION PHIBBY1,PHIBBY2,PHIBBY3,PHIBBY4,PHIBBY5,PHIBBY6,
     &  PHIBBY7
      DOUBLE PRECISION DUM

      INTEGER I

C     DATA CLIGHT/2.99792458D8/
C     DATA EMASSGEV/511003.3732832001D-9/

C     DIMENSION XLIM(12)

      Y=Y
      Z=Z

      X=XX
CERR18.3.93 IF (IBSYM.NE.0) X=DABS(XX)

      IF (IBFELD.NE.1) THEN   !FIRST CALL

       IF(IKBFORM.NE.0) THEN
          CALL BFORM
       ENDIF

       XLIM(1) =XM1-YSOFT1(1)
       XLIM(2) =XP1+YSOFT1(2)
       XLIM(3) =XM2-YSOFT2(1)
       XLIM(4) =XP2+YSOFT2(2)
       XLIM(5) =XM3-YSOFT3(1)
       XLIM(6) =XP3+YSOFT3(2)
       XLIM(7) =XM4-YSOFT4(1)
       XLIM(8) =XP4+YSOFT4(2)
       XLIM(9) =XM5-YSOFT5(1)
       XLIM(10)=XP5+YSOFT5(2)
       XLIM(11)=XM6-YSOFT6(1)
       XLIM(12)=XP6+YSOFT6(2)
       XLIM(13)=XM7-YSOFT7(1)
       XLIM(14)=XP7+YSOFT7(2)

       DO I=1,6
         IF(XLIM(I).GT.XLIM(I+1)) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)' *** ERROR IN BFELD ***'
            WRITE(LUNGFO,*)'MAGNETS COLLIDE, CHECK INPUT'
            WRITE(LUNGFO,*)'(NAMELIST BBFELD)'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN BFELD ***'
            WRITE(6,*)'MAGNETS COLLIDE, CHECK INPUT'
            WRITE(6,*)'(NAMELIST BBFELD)'
            WRITE(6,*)
            STOP
         ENDIF
       ENDDO
          IF (IBGAUSS.EQ.0) THEN

       XLBBY1=XP1-XM1
       XLBBY2=XP2-XM2
       XLBBY3=XP3-XM3
       XLBBY4=XP4-XM4
       XLBBY5=XP5-XM5
       XLBBY6=XP6-XM6
       XLBBY7=XP7-XM7

       PELEV=DSQRT(  (DMYENERGY-EMASSG1)*(DMYENERGY+EMASSG1)  )*1.D9

       IF (BBY1.NE.0.0) R0BBY1=PELEV/(CLIGHT1*BBY1)
       IF (BBY2.NE.0.0) R0BBY2=PELEV/(CLIGHT1*BBY2)
       IF (BBY3.NE.0.0) R0BBY3=PELEV/(CLIGHT1*BBY3)
       IF (BBY4.NE.0.0) R0BBY4=PELEV/(CLIGHT1*BBY4)
       IF (BBY5.NE.0.0) R0BBY5=PELEV/(CLIGHT1*BBY5)
       IF (BBY6.NE.0.0) R0BBY6=PELEV/(CLIGHT1*BBY6)
       IF (BBY7.NE.0.0) R0BBY7=PELEV/(CLIGHT1*BBY7)

       IF (R0BBY1.NE.0) THEN
          DUM=XLBBY1/R0BBY1
          IF (DABS(DUM).LE.1.D0) THEN
         PHIBBY1=DASIN(DUM)
          ELSE
         PHIBBY1=-9999.D0
          ENDIF
       ENDIF
       IF (R0BBY2.NE.0) THEN
          DUM=XLBBY2/R0BBY2
          IF (DABS(DUM).LE.1.D0) THEN
         PHIBBY2=DASIN(DUM)
          ELSE
         PHIBBY2=-9999.D0
          ENDIF
       ENDIF
       IF (R0BBY3.NE.0) THEN
          DUM=XLBBY3/R0BBY3
          IF (DABS(DUM).LE.1.D0) THEN
         PHIBBY3=DASIN(DUM)
          ELSE
         PHIBBY3=-9999.D0
          ENDIF
       ENDIF
       IF (R0BBY4.NE.0) THEN
          DUM=XLBBY4/R0BBY4
          IF (DABS(DUM).LE.1.D0) THEN
         PHIBBY4=DASIN(DUM)
          ELSE
         PHIBBY4=-9999.D0
          ENDIF
       ENDIF
       IF (R0BBY5.NE.0) THEN
          DUM=XLBBY5/R0BBY5
          IF (DABS(DUM).LE.1.D0) THEN
         PHIBBY5=DASIN(DUM)
          ELSE
         PHIBBY5=-9999.D0
          ENDIF
       ENDIF
       IF (R0BBY6.NE.0) THEN
          DUM=XLBBY6/R0BBY6
          IF (DABS(DUM).LE.1.D0) THEN
         PHIBBY6=DASIN(DUM)
          ELSE
         PHIBBY6=-9999.D0
          ENDIF
       ENDIF
       IF (R0BBY7.NE.0) THEN
         DUM=XLBBY7/R0BBY7
         IF (DABS(DUM).LE.1.D0) THEN
           PHIBBY7=DASIN(DUM)
         ELSE
           PHIBBY7=-9999.D0
         ENDIF
       ENDIF

       WRITE(LUNGFO,*)
       WRITE(LUNGFO,*)'      SR BFELD:'
       IF(IKBFORM.NE.0) THEN
          CALL BFORM
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'B0FORM,B0LP:      ',SNGL(B0FORM),SNGL(B0LP)
          WRITE(LUNGFO,*)'FB0N,FB0M  :      ',SNGL(FB0N),SNGL(FB0M)
          WRITE(LUNGFO,*)
       ENDIF
       WRITE(LUNGFO,*)
       WRITE(LUNGFO,*)'      Input parameters of XMi, XPi, YSOFTi, length XPi-XMi, magnetic field, bending radius, bending angle'
       WRITE(LUNGFO,*)'      (fringe fields not taken into account):'
       WRITE(LUNGFO,*)
       WRITE(LUNGFO,*)'      XM1,XP1,YSOFT1(1:2):',
     &                         SNGL(XM1),SNGL(XP1)
     &                        ,SNGL(YSOFT1(1)),SNGL(YSOFT1(2))
       WRITE(LUNGFO,*)'      XLBBY1,BBY1,R0BBY1,PHIBBY1:',
     &                   SNGL(XLBBY1),SNGL(BBY1),SNGL(R0BBY1),SNGL(PHIBBY1)
       WRITE(LUNGFO,*)'      XM2,XP2,YSOFT2(1:2):',
     &                         SNGL(XM2),SNGL(XP2)
     &                        ,SNGL(YSOFT1(1)),SNGL(YSOFT2(2))
       WRITE(LUNGFO,*)'      XLBBY2,BBY2,R0BBY2,PHIBBY2:',
     &                   SNGL(XLBBY2),SNGL(BBY2),SNGL(R0BBY2),SNGL(PHIBBY2)
       WRITE(LUNGFO,*)'      XM3,XP3,YSOFT3(1:2):',
     &                         SNGL(XM3),SNGL(XP3)
     &                        ,SNGL(YSOFT3(1)),SNGL(YSOFT3(2))
       WRITE(LUNGFO,*)'      XLBBY3,BBY3,R0BBY3,PHIBBY3:',
     &                   SNGL(XLBBY3),SNGL(BBY3),SNGL(R0BBY3),SNGL(PHIBBY3)
       WRITE(LUNGFO,*)'      XM4,XP4,YSOFT4(1:2):',
     &                         SNGL(XM4),SNGL(XP4)
     &                        ,SNGL(YSOFT4(1)),SNGL(YSOFT4(2))
       WRITE(LUNGFO,*)'      XLBBY4,BBY4,R0BBY4,PHIBBY4:',
     &                   SNGL(XLBBY4),SNGL(BBY4),SNGL(R0BBY4),SNGL(PHIBBY4)
       WRITE(LUNGFO,*)'      XM5,XP5,YSOFT5(1:2):',
     &                         SNGL(XM5),SNGL(XP5)
     &                        ,SNGL(YSOFT5(1)),SNGL(YSOFT5(2))
       WRITE(LUNGFO,*)'      XLBBY5,BBY5,R0BBY5,PHIBBY5:',
     &                   SNGL(XLBBY5),SNGL(BBY5),SNGL(R0BBY5),SNGL(PHIBBY5)
       WRITE(LUNGFO,*)'      XM6,XP6,YSOFT6(1:2):',
     &                         SNGL(XM6),SNGL(XP6)
     &                        ,SNGL(YSOFT6(1)),SNGL(YSOFT6(2))
       WRITE(LUNGFO,*)'      XLBBY6,BBY6,R0BBY6,PHIBBY6:',
     &                   SNGL(XLBBY6),SNGL(BBY6),SNGL(R0BBY6),SNGL(PHIBBY6)
       WRITE(LUNGFO,*)'      XM7,XP7,YSOFT7(1:2):',
     &   SNGL(XM7),SNGL(XP7)
     &   ,SNGL(YSOFT7(1)),SNGL(YSOFT7(2))
       WRITE(LUNGFO,*)'      XLBBY7,BBY7,R0BBY7,PHIBBY7:',
     &   SNGL(XLBBY7),SNGL(BBY7),SNGL(R0BBY7),SNGL(PHIBBY7)
      ENDIF

          IBFELD=1

      ENDIF

      BX=0.D0
      BY=0.D0
      BZ=0.D0

C=================================================== 17.10.90
      IF(IBGAUSS.NE.0) GOTO 1000
C=================================================== 17.10.90
CC    slam=(XP3-XM3)/(3.141592654D0*2.d0)
CC    BY=BBY3*DCOS((X-0.5D0*(XM3+XP3))/slam)*dcosh(y/slam)
CC    bx=-bby3*dsinh(y/slam)*dsin((x-0.5d0*(xm3+xp3))/slam)
CC    goto 999
C------------------------------
      IF (
     &  X.LT.(XM1-YSOFT1(1))
     &  .OR.
     &  (X.GT.(XP1+YSOFT1(2)) .AND. X.LT.(XM2-YSOFT2(1)))
     &  .OR.
     &  (X.GT.(XP2+YSOFT2(2)) .AND. X.LT.(XM3-YSOFT3(1)))
     &  .OR.
     &  (X.GT.(XP3+YSOFT3(2)) .AND. X.LT.(XM4-YSOFT4(1)))
     &  .OR.
     &  (X.GT.(XP4+YSOFT4(2)) .AND. X.LT.(XM5-YSOFT5(1)))
     &  .OR.
     &  (X.GT.(XP5+YSOFT5(2)) .AND. X.LT.(XM6-YSOFT6(1)))
     &  .OR.
     &  (X.GT.(XP6+YSOFT6(2)) .AND. X.LT.(XM7-YSOFT7(1)))
     &  .OR.
     &  (X.GT.(XP7+YSOFT7(2)))) goto 999

C--- 1. MAGNET

      IF (X.LE.XM1) THEN
       IF (YSOFT1(1).LT.1.D-31) goto 999
       BY=BBY1*0.5D0*(1.D0+DCOSD((X-XM1)/YSOFT1(1)*180.D0))
       goto 999
      ENDIF

      IF (X.LE.XP1) THEN
       BY=BBY1
       goto 999
      ENDIF

      IF (X.LE.(XP1+YSOFT1(2))) THEN
       IF (YSOFT1(2).LT.1.D-31) goto 999
       BY=BBY1*0.5D0*(1.D0+DCOSD((X-XP1)/YSOFT1(2)*180.D0))
       goto 999
      ENDIF

C--- 2. MAGNET

      IF (X.LE.XM2) THEN
       IF (YSOFT2(1).LT.1.D-31) goto 999
       BY=BBY2*0.5D0*(1.D0+DCOSD((X-XM2)/YSOFT2(1)*180.D0))
       goto 999
      ENDIF

      IF (X.LE.XP2) THEN
       BY=BBY2
       goto 999
      ENDIF

      IF (X.LE.(XP2+YSOFT2(2))) THEN
       IF (YSOFT2(2).LT.1.D-31) goto 999
       BY=BBY2*0.5D0*(1.D0+DCOSD((X-XP2)/YSOFT2(2)*180.D0))
       goto 999
      ENDIF

C--- 3. MAGNET

      IF (X.LE.XM3) THEN
       IF (YSOFT3(1).LT.1.D-31) goto 999
       BY=BBY3*0.5D0*(1.D0+DCOSD((X-XM3)/YSOFT3(1)*180.D0))
       goto 999
      ENDIF

      IF (X.LE.XP3) THEN
       BY=BBY3
       goto 999
      ENDIF

      IF (X.LE.(XP3+YSOFT3(2))) THEN
       IF (YSOFT3(2).LT.1.D-31) goto 999
       BY=BBY3*0.5D0*(1.D0+DCOSD((X-XP3)/YSOFT3(2)*180.D0))
       goto 999
      ENDIF

C--- 4. MAGNET

      IF (X.LE.XM4) THEN
       IF (YSOFT4(1).LT.1.D-31) goto 999
       BY=BBY4*0.5D0*(1.D0+DCOSD((X-XM4)/YSOFT4(1)*180.D0))
       goto 999
      ENDIF

      IF (X.LE.XP4) THEN
       BY=BBY4
       goto 999
      ENDIF

      IF (X.LE.(XP4+YSOFT4(2))) THEN
       IF (YSOFT4(2).LT.1.D-31) goto 999
       BY=BBY4*0.5D0*(1.D0+DCOSD((X-XP4)/YSOFT4(2)*180.D0))
       goto 999
      ENDIF

C--- 5. MAGNET

      IF (X.LE.XM5) THEN
       IF (YSOFT5(1).LT.1.D-31) goto 999
       BY=BBY5*0.5D0*(1.D0+DCOSD((X-XM5)/YSOFT5(1)*180.D0))
       goto 999
      ENDIF

      IF (X.LE.XP5) THEN
       BY=BBY5
       goto 999
      ENDIF

      IF (X.LE.(XP5+YSOFT5(2))) THEN
       IF (YSOFT5(2).LT.1.D-31) goto 999
       BY=BBY5*0.5D0*(1.D0+DCOSD((X-XP5)/YSOFT5(2)*180.D0))
       goto 999
      ENDIF

C--- 6. MAGNET

      IF (X.LE.XM6) THEN
        IF (YSOFT6(1).LT.1.D-31) goto 999
        BY=BBY6*0.5D0*(1.D0+DCOSD((X-XM6)/YSOFT6(1)*180.D0))
       goto 999
      ENDIF

      IF (X.LE.XP6) THEN
        BY=BBY6
       goto 999
      ENDIF

      IF (X.LE.(XP6+YSOFT6(2))) THEN
        IF (YSOFT6(2).LT.1.D-31) goto 999
        BY=BBY6*0.5D0*(1.D0+DCOSD((X-XP6)/YSOFT6(2)*180.D0))
       goto 999
      ENDIF

C--- 7. MAGNET

      IF (X.LE.XM7) THEN
        IF (YSOFT7(1).LT.1.D-31) goto 999
        BY=BBY7*0.5D0*(1.D0+DCOSD((X-XM7)/YSOFT7(1)*180.D0))
       goto 999
      ENDIF

      IF (X.LE.XP7) THEN
        BY=BBY7
       goto 999
      ENDIF

      IF (X.LE.(XP7+YSOFT7(2))) THEN
        IF (YSOFT7(2).LT.1.D-31) goto 999
        BY=BBY7*0.5D0*(1.D0+DCOSD((X-XP7)/YSOFT7(2)*180.D0))
       goto 999
      ENDIF


      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'*** ERROR IN BFELD ***'
      WRITE(LUNGFO,*)'SOMETHING IS WRONG, CHECK INPUT'
      WRITE(LUNGFO,*)
      WRITE(6,*)
      WRITE(6,*)'*** ERROR IN BFELD ***'
      WRITE(6,*)'SOMETHINGS WRONG, CHECK INPUT'
      WRITE(6,*)
         STOP
C=================================================== 17.10.90
C     !BERECHNE MAGNETFELD ALS UEBERLAGERUNG VON 3 GAUSSIANS
C     ! HSR-88-2, SEP.14, 1988, ODER TSKUBA KONF 1988 S.1851
1000    CONTINUE
      BY=0
      DO I=1,7
      BY=BY+
     &  B0GAUSS(I)*DEXP (
     &  -0.5D0*( (X-X0GAUSS(I)) / BSGGAUS(I) )**2)

C151091     &  -2.77259D0*( (X-X0GAUSS(I)) / BSGGAUS(I) )**2)
CERROR 18.10.90 SIEH LOGBUCH "PROBLEME"   BY=BY+
CERROR     &  B0GAUSS(I)*DEXP (
CERROR     &  -2.77259D2* (X-X0GAUSS(I))**2 / BSGGAUS(I) )
      ENDDO

C=================================================== 17.10.90

CERR18.3.93 999   IF (IBSYM.NE.0.AND.XX.LT.0.0) BX=-BX   !12.6.92
999   CONTINUE    ! BX=-BX IN MYBFELD
      RETURN


      END
