*CMZ :  4.00/14 22/12/2021  16.52.22  by  Michael Scheer
*CMZ :  4.00/11 19/05/2021  09.19.51  by  Michael Scheer
*CMZ :  3.05/10 08/08/2018  14.47.17  by  Michael Scheer
*CMZ :  3.02/03 23/10/2014  13.43.13  by  Michael Scheer
*CMZ :  3.01/00 06/05/2013  13.15.56  by  Michael Scheer
*CMZ :  3.00/01 19/03/2013  17.16.39  by  Michael Scheer
*CMZ :  3.00/00 14/03/2013  10.32.05  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  10.28.51  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.68/02 14/06/2012  13.07.05  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  09.55.58  by  Michael Scheer
*CMZ :  2.66/20 06/07/2011  10.03.53  by  Michael Scheer
*CMZ :  2.66/13 21/06/2010  15.29.34  by  Michael Scheer
*CMZ :  2.63/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.61/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.58/00 15/01/2007  09.59.51  by  Michael Scheer
*CMZ :  2.56/00 22/09/2005  15.19.02  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  10.48.09  by  Michael Scheer
*CMZ :  2.48/04 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.02  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  11.29.16  by  Michael Scheer
*CMZ :  2.20/05 13/03/2001  13.41.14  by  Michael Scheer
*CMZ :  2.16/08 31/10/2000  14.25.16  by  Michael Scheer
*CMZ :  2.16/06 29/08/2000  13.18.21  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.40.46  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  2.13/11 21/03/2000  13.34.41  by  Michael Scheer
*CMZ :  2.12/00 02/06/99  13.56.19  by  Michael Scheer
*CMZ :  1.03/06 25/06/98  17.04.16  by  Michael Scheer
*CMZ :  1.02/03 14/01/98  10.03.25  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.28  by  Michael Scheer
*CMZ : 00.01/12 11/09/96  17.54.03  by  Michael Scheer
*CMZ : 00.01/10 30/05/96  15.56.04  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.13.16  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.14  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.44  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WBTAB
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
*KEND.

C     WRITES MAGNETIC FIELD BY TO DATA FILE, FILE CAN BE READ FROM SR BTAB

      IMPLICIT NONE

*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,halbasy.
      include 'halbasy.cmn'
*KEND.

      INTEGER I,IT

      DOUBLE PRECISION EZ(3),EY(3),EZN,EYN
      DOUBLE PRECISION X,BX,BY,BZ,AX,AY,AZ,BI1,BI2,DX
      DOUBLE PRECISION BXU,BYU,BZU,BXL,BYL,BZL
      DOUBLE PRECISION AXM,AYM,AZM,AYINT,AZINT,AXP,AYP,AZP,BXM,BZM,BXP,BZP
      DOUBLE PRECISION BIQ,BYM,BYP
      DOUBLE PRECISION QMAX,SMAX
      DOUBLE PRECISION D2BZDZ,D2BYDZ,D2BZDY,D2BYDY
     &  ,DBYDZ,DBZDZ,VXZNNN
     &  ,DBYPDZ,DBYMDZ,DBZPDZ,DBZMDZ
     &  ,DBYUDY,DBYLDY
     &  ,DBZUDY,DBZLDY
     &  ,BEY,BEZ
     &  ,D52,D32,D21,D24
      DOUBLE PRECISION    DBY1DZ,D2BY1DZ
      DOUBLE PRECISION    DBZ1DZ,D2BZ1DZ
      DOUBLE PRECISION    D2BY1DY
      DOUBLE PRECISION    D2BZ1DY
      DOUBLE PRECISION    BY1,BZ1,BI11
      DOUBLE PRECISION BYDZS,BZDZS,BYDYS,BZDYS
      DOUBLE PRECISION BYDZS1,BZDZS1,BYDYS1,BZDYS1
      DOUBLE PRECISION BBYDZS,BBZDZS,BBYDYS,BBZDYS
      DOUBLE PRECISION VXN,VYN,VZN,V,DS
     &  ,DGAMMA

      DOUBLE PRECISION X2B(5),Y2B(5),Z2B(5)
     &  ,BFX(5),BFY(5),BFZ(5)
     &  ,AX2,AY2,AZ2
      DOUBLE PRECISION X1(5),Y1(5),Z1(5)
     &  ,VX1(5),VY1(5),VZ1(5),DT2,DT
     &  ,X2(5),Y2(5),Z2(5),VX2(5),VY2(5),VZ2(5)
     &  ,VXP(5),VYP(5),VZP(5)

c+self,IF=WBTABPARABEL.
      INTEGER IFAIL
      real (kind=16) XOPT,AOPT,XPAR(3),YPAR(3),APAR(3),YPPAR(3)
c+self.

c+self,IF=WBTABTUP.
      INTEGER NTUP_P,ICYCLE
      PARAMETER (NTUP_P=17)
      REAL*8 TUP_D(NTUP_P)
      CHARACTER(5) CHTAGS_D(NTUP_P)
      CHARACTER(6) CHTAGS(7)

      data chtags_d/'x','by','y','z','bi1','bi2'
     &  ,'bq','byzz','bzzz','byyy','bzyy','b2yzz'
     &  ,'b2zzz','b2yyy','b2zyy','ay','az'/

      data chtags/'x','y','z','d2bydy','d2bydz','d2bzdy','d2bzdz'/
c+self.

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,wbtab.
      include 'wbtab.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

c+self,IF=-WBTABTUP.
      REAL*4 SCALE
      DATA SCALE/1./
c+self.

      IF(BTABS.EQ.9999.) BTABS=XSTART
      IF(BTABE.EQ.9999.) BTABE=XSTOP
      IF(NPWBTAB.EQ.9999) NPWBTAB=NCO

      CALL hbookm(900,'WBTAB$',NTUP_P,'//WAVE',npwbtab,CHTAGS_D)
      CALL hbookm(901,'WBSEX$',7,'//WAVE',npwbtab,CHTAGS)

c+self,IF=-WBTABPARABEL.
C     WRITE(6,*)'*** WBTAB: SELECTION IN CMZ -WBTABPARABEL ***'
c+self.
      IF (FILEWBT.EQ.FILETB.and.irbtab.ne.0.or.irbtabzy.ne.0.
     &    or.irbtabxyz.ne.0.or.ifourbtabzy.ne.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN WBTAB ****'
        WRITE(LUNGFO,*)
     &    ' FILESNAMES FOR READING AND WRITING MAGNETIC FIELD DATA IDENTICAL '
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN WBTAB ****'
        WRITE(6,*)
     &    ' FILESNAMES FOR READING AND WRITING MAGNETIC FIELD DATA IDENTICAL '
        WRITE(6,*)
        STOP
      ENDIF

c+self,IF=-WBTABTUP
      OPEN(UNIT=LUNWBT,FILE=FILEWBT,STATUS='unknown')
      OPEN(UNIT=99,FILE='wbtab_sextupole.dat',STATUS='unknown')

      WRITE(LUNWBT,*) ICODE,CODE
      WRITE(LUNWBT,*) SCALE,SCALE

      WRITE(99,*) ICODE,CODE

      IF(IWBTAB.NE.100) THEN

        WRITE(LUNWBT,*) NPWBTAB

        WRITE(99,*) NPWBTAB
c+self,IF=WBTABTUP.
c      IF(IWBTAB.NE.100) THEN
c+self.
        BI1=0.0
        BI2=0.0
        BIQ=0.0
        BYDZS=0.0
        BZDZS=0.0
        BYDYS=0.0
        BZDYS=0.0
        BBYDZS=0.0
        BBZDZS=0.0
        BBYDYS=0.0
        BBZDYS=0.0
        QMAX=0.0
        SMAX=0.0
        AZINT=0.D0
        AYINT=0.D0

        X=BTABS
        DX=(BTABE-BTABS)/(NPWBTAB-1)

        DO I=1,NPWBTAB

          CALL MYBFELD(X,BTABY,BTABZ-BTABEPS,BXM,BYM,BZM,AXM,AYM,AZM)
          CALL MYBFELD(X,BTABY,BTABZ        ,BX ,BY ,BZ ,AX ,AY ,AZ )
          CALL MYBFELD(X,BTABY,BTABZ+BTABEPS,BXP,BYP,BZP,AXP,AYP,AZP)

          DBYPDZ=(BYP-BY)/BTABEPS
          DBYMDZ=(BY-BYM)/BTABEPS
          DBYDZ=(DBYPDZ+DBYMDZ)/2.
          D2BYDZ=(DBYPDZ-DBYMDZ)/BTABEPS

          DBZPDZ=(BZP-BZ)/BTABEPS
          DBZMDZ=(BZ-BZM)/BTABEPS
          DBZDZ=(DBZPDZ+DBZMDZ)/2.
          D2BZDZ=(DBZPDZ-DBZMDZ)/BTABEPS
          CALL MYBFELD(X,BTABY-BTABEPS,BTABZ,BXL,BYL,BZL,AX,AY,AZ)
          CALL MYBFELD(X,BTABY+BTABEPS,BTABZ,BXU,BYU,BZU,AX,AY,AZ)

          DBYUDY=(BYU-BY)/BTABEPS
          DBYLDY=(BY-BYL)/BTABEPS
          D2BYDY=(DBYUDY-DBYLDY)/BTABEPS

          DBZUDY=(BZU-BZ)/BTABEPS
          DBZLDY=(BZ-BZL)/BTABEPS
          D2BZDY=(DBZUDY-DBZLDY)/BTABEPS

          IF (I.EQ.1) THEN
            BY1=BY
            BZ1=BZ
            DBY1DZ=DBYDZ
            D2BY1DZ=D2BYDZ
            DBZ1DZ=DBZDZ
            D2BZ1DZ=D2BZDZ
            D2BY1DY=D2BYDY
            D2BZ1DY=D2BZDY
            BI11=0.D0
            BYDZS1=0.D0
            BZDZS1=0.D0
            BYDYS1=0.D0
            BZDYS1=0.D0
          ELSE IF (I.GT.1) THEN
            BI1=BI1+(BY1+BY)/2.D0*DX
            BI2=BI2+(BI11+BI1)/2.D0*DX
            BIQ=BIQ+(DBY1DZ+DBYDZ)/2.D0*DX

            BYDZS=BYDZS+(D2BY1DZ+D2BYDZ)/2.D0*DX
            BZDZS=BZDZS+(D2BZ1DZ+D2BZDZ)/2.D0*DX
            BYDYS=BYDYS+(D2BY1DY+D2BYDY)/2.D0*DX
            BZDYS=BZDYS+(D2BZ1DY+D2BZDY)/2.D0*DX

            BBYDZS=BBYDZS+(BYDZS1+BYDZS)/2.D0*DX
            BBZDZS=BBZDZS+(BZDZS1+BZDZS)/2.D0*DX
            BBYDYS=BBYDYS+(BYDYS1+BYDYS)/2.D0*DX
            BBZDYS=BBZDYS+(BZDYS1+BZDYS)/2.D0*DX

            AYINT=AYINT+(BZ1+BZ)/2.D0*DX
            AZINT=AZINT-(BY1+BY)/2.D0*DX
            BY1=BY
            BZ1=BZ
            DBY1DZ=DBYDZ
            BI11=BI1
            D2BY1DZ=D2BYDZ
            D2BZ1DZ=D2BZDZ
            D2BY1DY=D2BYDY
            D2BZ1DY=D2BZDY

            BYDZS1=BYDZS
            BZDZS1=BZDZS
            BYDYS1=BYDYS
            BZDYS1=BZDYS
          ENDIF

          IF(DABS(DBYDZ).GT.QMAX) QMAX=DBYDZ
          IF(DABS(D2BYDZ).GT.SMAX) SMAX=D2BYDZ

c+self,IF=WBTABTUP.
          TUP_D(1)=X
          TUP_D(2)=BY
          TUP_D(3)=BTABY
          TUP_D(4)=BTABZ
          TUP_D(5)=BI1
          TUP_D(6)=BI2
          TUP_D(7)=BIQ
          TUP_D(8)=BYDZS
          TUP_D(9)=BZDZS
          TUP_D(10)=BYDYS
          TUP_D(11)=BZDYS
          TUP_D(12)=BBYDZS
          TUP_D(13)=BBZDZS
          TUP_D(14)=BBYDYS
          TUP_D(15)=BBZDYS
          TUP_D(16)=AYINT
          TUP_D(17)=AZINT
          CALL hfm(900,TUP_D)
          TUP_D(1)=X
          TUP_D(2)=BTABY
          TUP_D(3)=BTABZ
          TUP_D(4)=D2BYDY
          TUP_D(5)=D2BYDZ
          TUP_D(6)=D2BZDY
          TUP_D(7)=D2BZDZ
          CALL hfm(901,TUP_D)
c+self,IF=-WBTABTUP.
          WRITE(99,'(7(1PE15.6))')SNGL(X),SNGL(BTABY),SNGL(BTABZ)
     &      ,SNGL(D2BYDY),SNGL(D2BYDZ)
     &      ,SNGL(D2BZDY),SNGL(D2BZDZ)
          WRITE(LUNWBT,'(2(1PE14.6),13(1PE11.3),2(1PE13.5))')
     &      SNGL(X),SNGL(BY)
     &      ,SNGL(BTABY),SNGL(BTABZ)
     &      ,SNGL(BI1),SNGL(BI2)
     &      ,SNGL(BIQ)
     &      ,SNGL(BYDZS)
     &      ,SNGL(BZDZS)
     &      ,SNGL(BYDYS)
     &      ,SNGL(BZDYS)
     &      ,SNGL(BBYDZS)
     &      ,SNGL(BBZDZS)
     &      ,SNGL(BBYDYS)
     &      ,SNGL(BBZDYS)
     &      ,SNGL(AYINT),SNGL(AZINT)
c+self.

          X=X+DX

        ENDDO !I

      ELSE  ! IWBTAB=100 I.E. ALONG TRAJECTORY

        DT=(XSTOP-XSTART)/NPWBTAB/CLIGHT1
c+self,if=-wbtabntup
        WRITE(LUNWBT,*) NPWBTAB
        WRITE(99,*) NPWBTAB
c+self.
        BI1=0.0
        BI2=0.0
        BIQ=0.0
        BYDZS=0.0
        BZDZS=0.0
        BYDYS=0.0
        BZDYS=0.0
        BBYDZS=0.0
        BBZDZS=0.0
        BBYDYS=0.0
        BBZDYS=0.0
        QMAX=0.0
        SMAX=0.0
        AZINT=0.D0
        AYINT=0.D0

        DT2=DT/2.D0

        X1(2)=WTRA(1,1,1)
        Y1(2)=WTRA(2,1,1)
        Z1(2)=WTRA(3,1,1)

        VX1(2)=WTRA(1,2,1)
        VY1(2)=WTRA(2,2,1)
        VZ1(2)=WTRA(3,2,1)

        V=DSQRT(VX1(2)*VX1(2)+VY1(2)*VY1(2)+VZ1(2)*VZ1(2))
        VXN=VX1(2)/V
        VYN=VY1(2)/V
        VZN=VZ1(2)/V

        VXZNNN=DSQRT(VXN*VXN+VZN*VZN)

        X1(1)=X1(2)+BTABEPS*VZN/VXZNNN
        Y1(1)=Y1(2)
        Z1(1)=Z1(2)-BTABEPS*VXN/VXZNNN

        X1(3)=X1(2)-BTABEPS*VZN/VXZNNN
        Y1(3)=Y1(2)
        Z1(3)=Z1(2)+BTABEPS*VXN/VXZNNN

        VXZNNN=DSQRT(VXN*VXN+VYN*VYN)

        X1(5)=X1(2)-BTABEPS*VYN/VXZNNN
        Y1(5)=Y1(2)+BTABEPS*VXN/VXZNNN
        Z1(5)=Z1(2)

        X1(4)=X1(2)+BTABEPS*VYN/VXZNNN
        Y1(4)=Y1(2)-BTABEPS*VXN/VXZNNN
        Z1(4)=Z1(2)

        VX1(1)=VX1(2)
        VY1(1)=VY1(2)
        VZ1(1)=VZ1(2)
        VX1(3)=VX1(2)
        VY1(3)=VY1(2)
        VZ1(3)=VZ1(2)
        VX1(4)=VX1(2)
        VY1(4)=VY1(2)
        VZ1(4)=VZ1(2)
        VX1(5)=VX1(2)
        VY1(5)=VY1(2)
        VZ1(5)=VZ1(2)

        DO IT=1,5

          X2(IT)=X1(IT)
          Y2(IT)=Y1(IT)
          Z2(IT)=Z1(IT)

          VX2(IT)=VX1(IT)
          VY2(IT)=VY1(IT)
          VZ2(IT)=VZ1(IT)

        ENDDO  !IT

        DO I=1,NPWBTAB

          DO IT=1,5

            X1(IT)=X2(IT)
            Y1(IT)=Y2(IT)
            Z1(IT)=Z2(IT)

            VX1(IT)=VX2(IT)
            VY1(IT)=VY2(IT)
            VZ1(IT)=VZ2(IT)

            X2B(IT)=X1(IT)+VX1(IT)*DT2
            Y2B(IT)=Y1(IT)+VY1(IT)*DT2
            Z2B(IT)=Z1(IT)+VZ1(IT)*DT2

            CALL MYBFELD(X2B(IT),Y2B(IT),Z2B(IT)
     &        ,BFX(IT),BFY(IT),BFZ(IT),AX2,AY2,AZ2)
            CALL BMOVETAYL(X1(IT),Y1(IT),Z1(IT),VX1(IT),VY1(IT),VZ1(IT)
     &        ,BFX(IT),BFY(IT),BFZ(IT),DT,
     &        X2(IT),Y2(IT),Z2(IT),VX2(IT),VY2(IT),VZ2(IT)
     &        ,VXP(IT),VYP(IT),VZP(IT),DMYGAMMA,ICHARGE,BMOVECUT,IUSTEP,IENELOSS,DGAMMA)

          ENDDO   !IT

          IF(I.EQ.1.OR.I.EQ.NPWBTAB) THEN
            DS=DT*V/2.D0
          ELSE
            DS=DT*V
          ENDIF   !I

          EY(1)=X1(5)-X1(2)
          EY(2)=Y1(5)-Y1(2)
          EY(3)=Z1(5)-Z1(2)
          EYN=DSQRT(EY(1)*EY(1)+EY(2)*EY(2)+EY(3)*EY(3))
          EY(1)=EY(1)/EYN
          EY(2)=EY(2)/EYN
          EY(3)=EY(3)/EYN

          EZ(1)=X1(3)-X1(2)
          EZ(2)=Y1(3)-Y1(2)
          EZ(3)=Z1(3)-Z1(2)
          EZN=DSQRT(EZ(1)*EZ(1)+EZ(2)*EZ(2)+EZ(3)*EZ(3))
          EZ(1)=EZ(1)/EZN
          EZ(2)=EZ(2)/EZN
          EZ(3)=EZ(3)/EZN


          D32=DSQRT((X1(3)-X1(2))**2+(Y1(3)-Y1(2))**2+(Z1(3)-Z1(2))**2)
          D21=DSQRT((X1(1)-X1(2))**2+(Y1(1)-Y1(2))**2+(Z1(1)-Z1(2))**2)

          BEY=
     &      (BFX(3)-BFX(2))*EY(1)+
     &      (BFY(3)-BFY(2))*EY(2)+
     &      (BFZ(3)-BFZ(2))*EY(3)
          DBYPDZ=BEY/D32

          BEY=
     &      (BFX(2)-BFX(1))*EY(1)+
     &      (BFY(2)-BFY(1))*EY(2)+
     &      (BFZ(2)-BFZ(1))*EY(3)
          DBYMDZ=BEY/D21

          BEZ=
     &      (BFX(3)-BFX(2))*EZ(1)+
     &      (BFY(3)-BFY(2))*EZ(2)+
     &      (BFZ(3)-BFZ(2))*EZ(3)
          DBZPDZ=BEZ/D32

          BEZ=
     &      (BFX(2)-BFX(1))*EZ(1)+
     &      (BFY(2)-BFY(1))*EZ(2)+
     &      (BFZ(2)-BFZ(1))*EZ(3)
          DBZMDZ=BEZ/D21

          DBYDZ=(DBYPDZ+DBYMDZ)/2.
          D2BYDZ=(DBYPDZ-DBYMDZ)/((D32+D21)/2.D0)
          DBZDZ=(DBZPDZ+DBZMDZ)/2.
          D2BZDZ=(DBZPDZ-DBZMDZ)/((D32+D21)/2.D0)

          D52=DSQRT((X1(5)-X1(2))**2+(Y1(5)-Y1(2))**2+(Z1(5)-Z1(2))**2)
          D24=DSQRT((X1(4)-X1(2))**2+(Y1(4)-Y1(2))**2+(Z1(4)-Z1(2))**2)

          BEY=
     &      (BFX(5)-BFX(2))*EY(1)+
     &      (BFY(5)-BFY(2))*EY(2)+
     &      (BFZ(5)-BFZ(2))*EY(3)
          DBYUDY=BEY/D52

          BEY=
     &      (BFX(2)-BFX(4))*EY(1)+
     &      (BFY(2)-BFY(4))*EY(2)+
     &      (BFZ(2)-BFZ(4))*EY(3)
          DBYLDY=BEY/D24

          BEZ=
     &      (BFX(5)-BFX(2))*EZ(1)+
     &      (BFY(5)-BFY(2))*EZ(2)+
     &      (BFZ(5)-BFZ(2))*EZ(3)
          DBZUDY=BEZ/D52

          BEZ=
     &      (BFX(2)-BFX(4))*EZ(1)+
     &      (BFY(2)-BFY(4))*EZ(2)+
     &      (BFZ(2)-BFZ(4))*EZ(3)
          DBZLDY=BEZ/D24

          D2BYDY=(DBYUDY-DBYLDY)/((D52+D24)/2.D0)
          D2BZDY=(DBZUDY-DBZLDY)/((D52+D24)/2.D0)

          IF (I.EQ.1) THEN
            BY1=BFY(2)
            BZ1=BFZ(2)
            DBY1DZ=DBYDZ
            D2BY1DZ=D2BYDZ
            DBZ1DZ=DBZDZ
            D2BZ1DZ=D2BZDZ
            D2BY1DY=D2BYDY
            D2BZ1DY=D2BZDY
            BI11=0.D0
            BYDZS1=0.D0
            BZDZS1=0.D0
            BYDYS1=0.D0
            BZDYS1=0.D0
          ELSE IF (I.GT.1) THEN
            BI1=BI1+(BY1+BFY(2))/2.D0*DS
            BI2=BI2+(BI11+BI1)/2.D0*DS
            BIQ=BIQ+(DBY1DZ+DBYDZ)/2.D0*DS

            BYDZS=BYDZS+(D2BY1DZ+D2BYDZ)/2.D0*DS
            BZDZS=BZDZS+(D2BZ1DZ+D2BZDZ)/2.D0*DS
            BYDYS=BYDYS+(D2BY1DY+D2BYDY)/2.D0*DS
            BZDYS=BZDYS+(D2BZ1DY+D2BZDY)/2.D0*DS

            BBYDZS=BBYDZS+(BYDZS1+BYDZS)/2.D0*DS
            BBZDZS=BBZDZS+(BZDZS1+BZDZS)/2.D0*DS
            BBYDYS=BBYDYS+(BYDYS1+BYDYS)/2.D0*DS
            BBZDYS=BBZDYS+(BZDYS1+BZDYS)/2.D0*DS

            AYINT=AYINT+(BZ1+BFZ(2))/2.D0*DS
            AZINT=AZINT-(BY1+BFY(2))/2.D0*DS

            BY1=BFY(2)
            BZ1=BFZ(2)
            DBY1DZ=DBYDZ
            BI11=BI1
            D2BY1DZ=D2BYDZ
            D2BZ1DZ=D2BZDZ
            D2BY1DY=D2BYDY
            D2BZ1DY=D2BZDY

            BYDZS1=BYDZS
            BZDZS1=BZDZS
            BYDYS1=BYDYS
            BZDYS1=BZDYS
          ENDIF

          IF(DABS(DBYDZ).GT.QMAX) QMAX=DBYDZ
          IF(DABS(D2BYDZ).GT.SMAX) SMAX=D2BYDZ

c+self,IF=WBTABTUP.
          TUP_D(1)=X1(2)
          TUP_D(2)=BFY(2)
          TUP_D(3)=Y1(2)
          TUP_D(4)=Z1(2)
          TUP_D(5)=BI1
          TUP_D(6)=BI2
          TUP_D(7)=BIQ
          TUP_D(8)=BYDZS
          TUP_D(9)=BZDZS
          TUP_D(10)=BYDYS
          TUP_D(11)=BZDYS
          TUP_D(12)=BBYDZS
          TUP_D(13)=BBZDZS
          TUP_D(14)=BBYDYS
          TUP_D(15)=BBZDYS
          TUP_D(16)=AYINT
          TUP_D(17)=AZINT
          CALL hfm(900,TUP_D)
          TUP_D(1)=X1(2)
          TUP_D(2)=Y1(2)
          TUP_D(3)=Z1(2)
          TUP_D(4)=D2BYDY
          TUP_D(5)=D2BYDZ
          TUP_D(6)=D2BZDY
          TUP_D(7)=D2BZDZ
          CALL hfm(901,TUP_D)
c+self,IF=-WBTABTUP.
          WRITE(99,'(7(1PE15.6))')SNGL(X1(2)),SNGL(Y1(2)),SNGL(Z1(2))
     &      ,SNGL(D2BYDY),SNGL(D2BYDZ)
     &      ,SNGL(D2BZDY),SNGL(D2BZDZ)
          WRITE(LUNWBT,'(2(1PE14.6),13(1PE11.3),2(1PE13.5))')
     &      SNGL(X1(2)),SNGL(BFY(2))
     &      ,SNGL(Y1(2)),SNGL(Z1(2))
     &      ,SNGL(BI1),SNGL(BI2)
     &      ,SNGL(BIQ)
     &      ,SNGL(BYDZS)
     &      ,SNGL(BZDZS)
     &      ,SNGL(BYDYS)
     &      ,SNGL(BZDYS)
     &      ,SNGL(BBYDZS)
     &      ,SNGL(BBZDZS)
     &      ,SNGL(BBYDYS)
     &      ,SNGL(BBZDYS)
     &      ,SNGL(AYINT),SNGL(AZINT)

c+self.
        ENDDO  !NPWBTAB

      ENDIF   !IWBTAB

      CLOSE(LUNWBT)
      CLOSE(99)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR WBTAB: '
      WRITE(LUNGFO,*)
      IF (KHALBASY.NE.0.AND.IAHWFOUR.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN WBTAB:'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'FOR ASYMETIC HALBACH-MODEL SEXTUPOL-TERMS ARE WRONG'
        WRITE(LUNGFO,*)'SINCE IAHWFOUR=0'
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN WBTAB:'
        WRITE(6,*)
        WRITE(6,*)'FOR ASYMETIC HALBACH-MODEL SEXTUPOL-TERMS ARE WRONG'
        WRITE(6,*)'SINCE IAHWFOUR=0'
      ENDIF
      IF (IRFILF.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN WBTAB:'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    'FOR FOURIER-EXPANDED FIELDS SEXTUPOLE-TERMS ARE WRONG'
        WRITE(LUNGFO,*)'SINCE 3D-EFFECTS ARE NOT TAKEN INTO ACCOUNT, I.E.'
        WRITE(LUNGFO,*)'Kx=CONST.'
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN WBTAB:'
        WRITE(6,*)
        WRITE(6,*)
     &    'FOR FOURIER-EXPANDED FIELDS SEXTUPOLE-TERMS ARE WRONG'
        WRITE(6,*)'SINCE 3D-EFFECTS ARE NOT TAKEN INTO ACCOUNT, I.E.'
        WRITE(6,*)'Kx=CONST.'
      ENDIF
      WRITE(LUNGFO,*)
c+self,IF=-WBTABTUP.
      WRITE(LUNGFO,*)'     written to file:'
      WRITE(LUNGFO,*)'     ',FILEWBT
      WRITE(LUNGFO,*)
     &  '     For exact format please refer to wave.in, namelist WBTABN'
c+self,IF=WBTABTUP.
      WRITE(LUNGFO,*)'     written to NTUPLE'
      CALL MHROUT(900,ICYCLE,' ')
      CALL hdeletm(900)
      CALL MHROUT(901,ICYCLE,' ')
      CALL hdeletm(901)
c+self.
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     x-intervall:',SNGL(BTABS),SNGL(BTABE)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     first field integral (T m):                ',
     &  SNGL(BI1)
      WRITE(LUNGFO,*)'     second field integral (T m**2):            ',
     &  SNGL(BI2)
      WRITE(LUNGFO,*)'     quadrupole like integral (T):              ',
     &  SNGL(BIQ)
      WRITE(LUNGFO,*)'     maximum gradient (T/m):                    ',
     &  SNGL(QMAX)
      WRITE(LUNGFO,*)'     sextupole like integral d2By/dz**2 (T/m):  ',
     &  SNGL(BYDZS)
      WRITE(LUNGFO,*)'     sextupole like integral d2By/dy**2 (T/m):  ',
     &  SNGL(BYDYS)
      WRITE(LUNGFO,*)'     maximum gradient (T/m**2):                 ',
     &  SNGL(SMAX)

      IF(IWBTAB.EQ.100) THEN
        WRITE(LUNGFO,*)'     (along trajectory)'
      ELSE
        WRITE(LUNGFO,*)
     &    '     (along straight line, z =',
     &    SNGL(BTABZ),' y =',SNGL(BTABY),')'
      ENDIF !IWBTAB
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

      RETURN
      END
