*CMZ :  4.00/04 27/08/2019  11.49.27  by  Michael Scheer
*CMZ :  3.04/00 23/01/2018  12.41.22  by  Michael Scheer
*CMZ :  3.03/04 18/12/2017  12.11.39  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  10.38.36  by  Michael Scheer
*CMZ :  2.63/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.61/01 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.61/00 30/01/2007  19.54.23  by  Michael Scheer
*-- Author :    Michael Scheer   30/01/2007
      SUBROUTINE IDTRMSHGF(XIN,YIN,ZIN,VXI,VYI,VZI,XF,YF,ZF,VXF,VYF,VZF,
     &  XCO0,YCO0,ZCO0,VXCO0,VYCO0,VZCO0,XCOF,YCOF,ZCOF,VXCOF,VYCOF,VZCOF,
     &  coefile)
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

C--- TRACKING WITH POLYNOMIAL-EXPANSION OF GENERATING FUNCTION

C XIN,YIN,ZIN,XF,YF,ZF ARE COORDINATES IN LAB

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,optic.
      include 'optic.cmn'
*KEEP,genfun.
      include 'genfun.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DOUBLE PRECISION
     &  EWS(3),EWY(3),EWZ(3),EN,
     &  EWSF(3),EWYF(3),EWZF(3),YW,ZW,ZPW,YPW,SPW,
     &  YIN,ZIN,YF,ZF,YPF,ZPF,YACC,X0,Y0,Z0,
     &  XIN,VXI,VYI,VZI,V0,VXINN,VYINN,VZINN,VXF,VYF,VZF,
     &  XCO0,YCO0,ZCO0,VXCO0,VYCO0,VZCO0,XCOF,YCOF,ZCOF,VXCOF,VYCOF,VZCOF,
     &  BXF,BYF,BZF,AXF,AYF,AZF,SNENN2,PAX,PAY,DXX,DYY,DFFX,FFX,FFY,DFFY,
     &  BRHOABS,SNENN,sin

      REAL*8 AKOEFF(NORDNG,NORDNG,NORDNG,NORDNG)
      REAL*8 QXPOW(NORDNG+1),PXPOW(NORDNG+1)
      REAL*8 QYPOW(NORDNG+1),PYPOW(NORDNG+1)
      REAL*8 TRANSFM(4,4)
      REAL*8 A1100,A2000,A0200,A1000,A0100,A0011,A0020,A0002
      REAL*8 PXI,QXI,QYI,PYI,PX,PY,QX,QY,DX,DY,DZ,YWF,ZWF,XF,
     &  XL,YL,ZL,RNENN,pxi0,pyi0,pxf0,pyf0

      REAL*8
     &   XF0,YF0,ZF0
     &  ,ZP0,YP0
     &  ,BX0,BY0,BZ0
     &  ,AX0,AY0,AZ0
     &  ,ZPF0,YPF0
     &  ,BXF0,BYF0,BZF0
     &  ,AXF0,AYF0,AZF0

      REAL*8 AXR0,AYR0,AZR0
      REAL*8 AXR,AYR,AZR
      REAL*8 BXI,BYI,BZI
      REAL*8 AXI,AYI,AZI
      REAL*8 AXRF0,AYRF0,AZRF0,A0SCALE,PEL

      INTEGER ICAL,LUNCOE,JCHARGE,MORDNG,JMAX,IREAD,JWRITE,IWRITE,
     &  I,J,K,L,JLOOP,IZAPER,IYAPER,
     &  I1,J1,K1,L1,
     &  I2,J2,K2,L2,IPOW,nkoef

      CHARACTER(*) COEFILE
      CHARACTER*65 CODETRA

      DATA LUNCOE/89/
      DATA ICAL/0/

      DATA QXPOW(1)/0./
      DATA PXPOW(1)/0./
      DATA QYPOW(1)/0./
      DATA PYPOW(1)/0./

      DATA QXPOW(2)/1./
      DATA PXPOW(2)/1./
      DATA QYPOW(2)/1./
      DATA PYPOW(2)/1./

      DATA JMAX/20/,YACC/1.D-10/

C--- INITIALIZATION {

      IF (ICAL.EQ.0) THEN

        PEL=EMASSE1*DSQRT( (DMYGAMMA+1.0D0)*(DMYGAMMA-1.0D0) )
        BRHOABS=PEL/CLIGHT1  !ABSOLUTE VALUE

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      IDTRMSHGF:'
        WRITE(LUNGFO,*)

c        COEFILE='wave_erzfun.in'

C Einlesen der Koeeficenten  GF(i,j) !!!

        OPEN(UNIT=LUNCOE,FILE = COEFILE,STATUS ='OLD')

        IREAD=0

        READ(LUNCOE,'(1A65)')CODETRA
        READ(LUNCOE,*)JCHARGE,MORDNG

        IF (JCHARGE.NE.ICHARGE) THEN
          WRITE(LUNGFO,*)'*** ERROR IN IDTRMSHGF:  ***'
          WRITE(LUNGFO,*)'GENERATING FUNCTION REFERES TO CHARGE',JCHARGE
          WRITE(6,*)'*** ERROR IN IDTRMSHGF:  ***'
          WRITE(6,*)'GENERATING FUNCTION REFERES TO CHARGE',JCHARGE
          STOP '*** PROGRAM WAVE ABORTED ***'
        ENDIF !JCHARGE

        IF (MORDNG.GT.NORDNG) THEN
          WRITE(LUNGFO,*)'*** ERROR IN IDTRMSHGF:  ***'
          WRITE(LUNGFO,*)'*** DIMENSION NORDNG EXCEEDED'
          WRITE(LUNGFO,*)
     &    'TOO HIGH! INCREASE PARAMETER NORDNG IN GENFUN.CMN'
          WRITE(6,*)'*** ERROR IN IDTRMSHGF:  ***'
          WRITE(6,*)'*** DIMENSION NORDNG EXCEEDED'
          WRITE(6,*)
     &    'TOO HIGH! INCREASE PARAMETER NORDNG IN GENFUN.CMN'
          STOP '*** PROGRAM WAVE ABORTED ***'
        ENDIF

        READ(LUNCOE,*)ZAPERT,YAPERT,DLAPER,sin

        READ(LUNCOE,*)X0,Y0,Z0
     &    ,ZP0,YP0
     &    ,BX0,BY0,BZ0
     &    ,AX0,AY0,AZ0

        READ(LUNCOE,*)XF0,YF0,ZF0
     &    ,ZPF0,YPF0
     &    ,BXF0,BYF0,BZF0
     &    ,AXF0,AYF0,AZF0

        READ(LUNCOE,*)EWS(1),EWS(2),
     &    EWS(3)
        READ(LUNCOE,*)EWSF(1),EWSF(2),
     &    EWSF(3)
        READ(LUNCOE,*)A0SCALE

        IF (
     &      ABS(X0-XCO0).GT.YACC.OR.
     &      ABS(Y0-YCO0).GT.YACC.OR.
     &      ABS(Z0-ZCO0).GT.YACC.OR.
     &      ABS(YP0-VYCO0/VXCO0).GT.YACC.OR.
     &      ABS(ZP0-VZCO0/VXCO0).GT.YACC.OR.
     &      ABS(XF0-XCOF).GT.YACC.OR.
     &      ABS(YF0-YCOF).GT.YACC.OR.
     &      ABS(ZF0-ZCOF).GT.YACC.OR.
     &      ABS(YPF0-VYCOF/VXCOF).GT.YACC.OR.
     &      ABS(ZPF0-VZCOF/VXCOF).GT.YACC
     &      ) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '*** WARNING IN IDTRMSHGF: DISCREPANCIES OF CLOSED ORBITS'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)' REFERENCE ORBIT AT ENTRANCE:'
          WRITE(LUNGFO,*)XCO0
          WRITE(LUNGFO,*)YCO0
          WRITE(LUNGFO,*)ZCO0
          WRITE(LUNGFO,*)VYCO0/VXCO0
          WRITE(LUNGFO,*)VZCO0/VXCO0
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)' REFERENCE ORBIT AT EXIT:'
          WRITE(LUNGFO,*)XCOF
          WRITE(LUNGFO,*)YCOF
          WRITE(LUNGFO,*)ZCOF
          WRITE(LUNGFO,*)VYCOF/VXCOF
          WRITE(LUNGFO,*)VZCOF/VXCOF
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)' CLOSED ORBIT OF GF AT ENTRANCE:'
          WRITE(LUNGFO,*)X0
          WRITE(LUNGFO,*)Y0
          WRITE(LUNGFO,*)Z0
          WRITE(LUNGFO,*)YP0
          WRITE(LUNGFO,*)ZP0
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)' CLOSED ORBIT OF GF AT EXIT:'
          WRITE(LUNGFO,*)XF0
          WRITE(LUNGFO,*)YF0
          WRITE(LUNGFO,*)ZF0
          WRITE(LUNGFO,*)YPF0
          WRITE(LUNGFO,*)ZPF0
        ENDIF

11      CONTINUE

        READ (LUNCOE,*,END=99) I,J,K,L
     &    ,AKOEFF(I+1,J+1,K+1,L+1)

        IF(I+1.GT.NORDNG.OR.J+1.GT.NORDNG
     &      .OR.K+1.GT.NORDNG.OR.L+1.GT.NORDNG) THEN
          WRITE(6,*)
          WRITE(6,*)'*** ERROR SR IDTRMSHGF ***'
          WRITE(6,*)'ORDER OF COEFFICIENTS ON FILE'
          WRITE(6,*)coefile(1:len_trim(coefile))
          WRITE(6,*)
     &      'TOO HIGH INCREASE PARAMETER NORDNG IN GENFUN.CMN'
          WRITE(6,*)
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR SR IDTRMSHGF ***'
          WRITE(LUNGFO,*)'ORDER OF COEFFICIENTS ON FILE'
          WRITE(LUNGFO,*)coefile(1:len_trim(coefile))
          WRITE(LUNGFO,*)
     &      'TOO HIGH INCREASE PARAMETER NORDNG IN GENFUN.CMN'
          WRITE(LUNGFO,*)
          STOP
        ENDIF
        IREAD=IREAD+1
        GOTO 11

99      CONTINUE

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      ',IREAD,' coeffs. read from file:'
        WRITE(LUNGFO,*)coefile(1:len_trim(coefile))
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)CODETRA
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)

        CLOSE(LUNCOE)

C--- CANONICAL MOMENTUM OF REFERENCE ORBIT

C     UNIT-VECTOR EWZ=[EWS,(0,1,0)] (CROSS-PRODUCT)

        EN=1.D0/DSQRT(
     &    EWS(3)*EWS(3)
     &    +EWS(1)*EWS(1))

        EWZ(1)=-EWS(3)*EN
        EWZ(2)= 0.
        EWZ(3)=+EWS(1)*EN

C     UNIT-VECTOR EWY=[EWZ,EWS]

        EWY(1)=-EWS(2)*EWZ(3)
     &    + EWS(3)*EWZ(2)
        EWY(2)=-EWS(3)*EWZ(1)
     &    + EWS(1)*EWZ(3)
        EWY(3)=-EWS(1)*EWZ(2)
     &    + EWS(2)*EWZ(1)

        EN=1.0D0/DSQRT(
     &    EWSF(3)*EWSF(3)
     &    +EWSF(1)*EWSF(1))
        EWZF(1)=-EWSF(3)*EN
        EWZF(2)= 0.
        EWZF(3)=+EWSF(1)*EN

        EWYF(1)=-EWSF(2)*EWZF(3)
     &    + EWSF(3)*EWZF(2)
        EWYF(2)=-EWSF(3)*EWZF(1)
     &    + EWSF(1)*EWZF(3)
        EWYF(3)=-EWSF(1)*EWZF(2)
     &    + EWSF(2)*EWZF(1)

        AXR0=AX0*EWS(1)
     &    +AY0*EWS(2)
     &    +AZ0*EWS(3)

        AYR0=AX0*EWY(1)
     &    +AY0*EWY(2)
     &    +AZ0*EWY(3)

        AZR0=AX0*EWZ(1)
     &    +AY0*EWZ(2)
     &    +AZ0*EWZ(3)

        AXRF0=AXF0*EWSF(1)
     &    +AYF0*EWSF(2)
     &    +AZF0*EWSF(3)

        AYRF0=AXF0*EWYF(1)
     &    +AYF0*EWYF(2)
     &    +AZF0*EWYF(3)

        AZRF0=AXF0*EWZF(1)
     &    +AYF0*EWZF(2)
     &    +AZF0*EWZF(3)

c18.1.2018{
        if (icharge.gt.0) then
          RNENN=SQRT(1.0D0 + ZP0**2 + YP0**2)
          PXI0=AZR0/BRHOABS + ZP0/RNENN
          PYI0=AYR0/BRHOABS + YP0/RNENN
          RNENN=SQRT(1.0D0 + ZPf0**2 + YPf0**2)
          PXf0=AZRf0/BRHOABS + ZPf0/RNENN
          PYf0=AYRf0/BRHOABS + YPf0/RNENN
        else
          RNENN=SQRT(1.0D0 + ZP0**2 + YP0**2)
          PXI0=-AZR0/BRHOABS + ZP0/RNENN
          PYI0=-AYR0/BRHOABS + YP0/RNENN
          RNENN=SQRT(1.0D0 + ZPf0**2 + YPf0**2)
          PXf0=-AZRf0/BRHOABS + ZPf0/RNENN
          PYf0=-AYRf0/BRHOABS + YPf0/RNENN
        endif
c18.1.2018}

C--- LINEARE TRANSFERMATRIX OHNE KOPPLUNG-TERME

        A1100=AKOEFF(2,2,1,1)
        A2000=AKOEFF(3,1,1,1)
        A0200=AKOEFF(1,3,1,1)
        A1000=AKOEFF(2,1,1,1)
        A0100=AKOEFF(1,2,1,1)

        A0011=AKOEFF(1,1,2,2)
        A0020=AKOEFF(1,1,3,1)
        A0002=AKOEFF(1,1,1,3)

        TRANSFM(1,1)=(-4.D0*AKOEFF(3,1,1,1)
     &    *AKOEFF(1,3,1,1)*AKOEFF(1,1,2,2)
     &    +2.D0*AKOEFF(3,1,1,1)*
     &    AKOEFF(1,2,2,1)*AKOEFF(1,2,1,2)
     &    +AKOEFF(2,2,1,1)**2.D0*AKOEFF(1,1,2,2)
     &    -AKOEFF(2,2,1,1)*
     &    AKOEFF(2,1,2,1)*AKOEFF(1,2,1,2)
     &    -AKOEFF(2,2,1,1)*AKOEFF(2,1,1,2)
     &    *AKOEFF(1,2,2,1)+2.D0*
     &    AKOEFF(2,1,2,1)*AKOEFF(2,1,1,2)*AKOEFF(1,3
     &    ,1,1))/(AKOEFF(2,2,1,1)
     &    *AKOEFF(1,1,2,2)-
     &    AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
        TRANSFM(1,2)=(2.D0*AKOEFF(1,3,1,1)
     &    *AKOEFF(1,1,2,2)-AKOEFF(1,2,2,1)
     &    *AKOEFF(1,2,1,2))/(AKOEFF
     &    (2,2,1,1)*AKOEFF(1,1,2,2)-AKOEFF(2,1,1,
     &    2)*AKOEFF(1,2,2,1))
        TRANSFM(1,3)=(AKOEFF(2,2,1,1)*AKOEFF(1,2,2,1
     &    )*AKOEFF(1,1,2,2)-2.D0
     &    *AKOEFF(2,2,1,1)*AKOEFF(
     &    1,2,1,2)*AKOEFF(1,1,3,1)
     &    -2.D0*AKOEFF(2,1,2,1)*AKOEFF(1,3,1,1)
     &    *AKOEFF(1,1,2,2)+AKOEFF(2,1,2,1)
     &    *AKOEFF(1,2,2,1)*AKOEFF(1,2,1,2)
     &    +4.D0*AKOEFF(2,1,1,2)*AKOEFF(1,3,1,1)*
     &    AKOEFF(1,1,3,1)-AKOEFF(2,1,1,2)
     &    *AKOEFF(1,2,2,1)**2)/(AKOEFF(2,2,1,1)
     &    *AKOEFF(1,1,2,2)-AKOEFF(2,1,1,2)
     &    *AKOEFF(1,2,2,1))
        TRANSFM(1,4)=(AKOEFF(2,2,1,1)*AKOEFF(1,2,1,2
     &    )-2.D0*AKOEFF(2,1,1,2)
     &    *AKOEFF(1,3,1,1))/(AKOEFF(2,2,1,1)
     &    *AKOEFF(1,1,2,2)-AKOEFF(2,1,1,
     &    2)*AKOEFF(1,2,2,1))
        TRANSFM(2,1)=(-2.D0*AKOEFF(3,1,1,1)
     &    *AKOEFF(1,1,2,2)+AKOEFF(2,1,2,1)
     &    *AKOEFF(2,1,1,2))/(
     &    AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)
     &    -AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
        TRANSFM(2,2)=AKOEFF(1,1,2,2)/
     &    (AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)
     &    -AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
        TRANSFM(2,3)=(-AKOEFF(2,1,2,1)*AKOEFF(1,1,2,
     &    2)+2.D0*AKOEFF(2,1,1,2)
     &    *AKOEFF(1,1,3,1))/(
     &    AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)
     &    -AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
        TRANSFM(2,4)=(-AKOEFF(2,1,1,2))/
     &    (AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)
     &    -AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
        TRANSFM(3,1)=(4.D0*AKOEFF(3,1,1,1)
     &    *AKOEFF(1,2,2,1)*AKOEFF(1,1,1,3)
     &    -2.D0*AKOEFF(3,1,1,1)*
     &    AKOEFF(1,2,1,2)*AKOEFF(1,1,2,2)
     &    -2.D0*AKOEFF(2,2,1,1)*AKOEFF(2,1,2,1)
     &    *AKOEFF(1,1,1,3)+
     &    AKOEFF(2,2,1,1)*AKOEFF(2,1,1,2)
     &    *AKOEFF(1,1,2,2)+AKOEFF(2,1,2,1)
     &    *AKOEFF(2,1,1,2)*AKOEFF(1,2,1,2)
     &    -AKOEFF(2,1,1,2)**2.D0*AKOEFF(1,
     &    2,2,1))/(AKOEFF(2,2,1,1)
     &    *AKOEFF(1,1,2,2)
     &    -AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
        TRANSFM(3,2)=(-2.D0*AKOEFF(1,2,2,1)
     &    *AKOEFF(1,1,1,3)+AKOEFF(1,2,1,2)
     &    *AKOEFF(1,1,2,2))/(AKOEFF(2,2,1,1)
     &    *AKOEFF(1,1,2,2)-AKOEFF(2,1
     &    ,1,2)*AKOEFF(1,2,2,1))
        TRANSFM(3,3)=(-4.D0*AKOEFF(2,2,1,1)
     &    *AKOEFF(1,1,
     &    3,1)*AKOEFF(1,1,1,3)
     &    +AKOEFF(2,2,1,1)*
     &    AKOEFF(1,1,2,2)**2+2.D0*AKOEFF(2,1,2,1)
     &    *AKOEFF(
     &    1,2,2,1)*AKOEFF(1,1,1,3)
     &    -AKOEFF(2,1,2,1)*AKOEFF(1,2,1,2)
     &    *AKOEFF(1,1,2,2)-AKOEFF(2,1,1,2)
     &    *AKOEFF(1,2,2,1)*AKOEFF(1,1,2,
     &    2)+2.D0*AKOEFF(2,1,1,2)
     &    *AKOEFF(1,2,1,2)*
     &    AKOEFF(1,1,3,1))/(AKOEFF(2,2,1,1)
     &    *AKOEFF(1,1,2,2)-AKOEFF(2,1,1,2)
     &    *AKOEFF(1,2,2,1))
        TRANSFM(3,4)=(2.D0*AKOEFF(2,2,1,1)
     &    *AKOEFF(1,1,1,3)-AKOEFF(2,1,1,2)
     &    *AKOEFF(1,2,1,2))/(AKOEFF
     &    (2,2,1,1)*AKOEFF(1,1,2,2)
     &    -AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
        TRANSFM(4,1)=(2.D0*AKOEFF(3,1,1,1)
     &    *AKOEFF(1,2,2,1)-AKOEFF(2,2,1,1)
     &    *AKOEFF(2,1,2,1))/(AKOEFF
     &    (2,2,1,1)*AKOEFF(1,1,2,2)
     &    -AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
        TRANSFM(4,2)=(-AKOEFF(1,2,2,1))
     &    /(AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)
     &    -AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
        TRANSFM(4,3)=(-2.D0*AKOEFF(2,2,1,1)
     &    *AKOEFF(1,1,3,1)+AKOEFF(2,1,2,1)
     &    *AKOEFF(1,2,2,1))/(
     &    AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)
     &    -AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
        TRANSFM(4,4)=AKOEFF(2,2,1,1)/
     &    (AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)
     &    -AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '      Linear Transfer Matrix from Generating Function'
        WRITE(LUNGFO,*)
        DO IWRITE=1,4
          WRITE(LUNGFO,*)'      ',(SNGL(TRANSFM(IWRITE,JWRITE))
     &      ,JWRITE=1,4)
        ENDDO
        WRITE(LUNGFO,*)

        ICAL=1

      ENDIF   !JCAL.EQ.0

C--- INITIALIZATION }

C --- ARE WE IN ENTRANCE PLAN?

      DX=XIN-X0
      DY=YIN-Y0
      DZ=ZIN-Z0

      IF (DX*EWS(1)+DY*EWS(2)+DZ*EWS(3) .GT. 1.0D-9) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '*** WARNING IN IDTRMSHGF: STARTING POINT NOT IN ENTRANCE PLANE'
        WRITE(LUNGFO,*)'XIN:',XIN
        WRITE(LUNGFO,*)'YIN:',YIN
        WRITE(LUNGFO,*)'ZIN:',ZIN
        WRITE(LUNGFO,*)
      ENDIF

      YW=DY
      ZW=DZ

      V0=SQRT(VXI**2+VYI**2+VZI**2)

      VXINN=VXI/V0
      VYINN=VYI/V0
      VZINN=VZI/V0

      SPW=VXINN*EWS(1)+VYINN*EWS(2)+VZINN*EWS(3)
      YPW=VXINN*EWY(1)+VYINN*EWY(2)+VZINN*EWY(3)
      ZPW=VXINN*EWZ(1)+VYINN*EWZ(2)+VZINN*EWZ(3)

      YPW=YPW/SPW
      ZPW=ZPW/SPW

      IF (IERZFUN.LT.0.0) THEN

C--- ONLY APPLY LINEAR TRANSFER MATRIX

        ZWF =   TRANSFM(1,1)*ZW + TRANSFM(1,2)*ZPW
     &    +   TRANSFM(1,3)*YW + TRANSFM(1,4)*YPW
        ZPF =  TRANSFM(2,1)*ZW + TRANSFM(2,2)*ZPW
     &    +   TRANSFM(2,3)*YW + TRANSFM(2,4)*YPW
        YWF =   TRANSFM(3,1)*ZW + TRANSFM(3,2)*ZPW
     &    +   TRANSFM(3,3)*YW + TRANSFM(3,4)*YPW
        YPF =  TRANSFM(4,1)*ZW + TRANSFM(4,2)*ZPW
     &    +   TRANSFM(4,3)*YW + TRANSFM(4,4)*YPW

        XF=XF0+ZWF*EWZF(1)+YWF*EWYF(1)
        YF=YF0+ZWF*EWZF(2)+YWF*EWYF(2)
        ZF=ZF0+ZWF*EWZF(3)+YWF*EWYF(3)

cmsh20171128 VXF=V0*EWSF(1)
cmsh20171128 VYF=VXF*YPF
cmsh20171128 VZF=VXF*ZPF

        VXF=V0*vxinn
        VYF=v0*vyinn+vxf*ypf
        VZF=v0*vzinn+vxf*zpf

        RETURN

      ENDIF !IERZFUN .LT.0

C--- HARDWARE APERTURE

      IF (DABS(ZW).GT.DABS(ZAPERT)
     &    .OR.  DABS(ZW+ZPW*DLAPER).GT.DABS(ZAPERT)) THEN
        IZAPER=1
        WRITE(6,*)'*** HORI. APERTURE WARNING IN IDTRMSHGF:'
        WRITE(6,*)'ZW,ZPW:',ZW,ZPW
      ENDIF

      IF (DABS(YW).GT.DABS(YAPERT)
     &    .OR.  DABS(YW+YPW*DLAPER).GT.DABS(YAPERT)) THEN
        IYAPER=1
        WRITE(6,*)'*** VERTICAL APERTURE WARNING IN IDTRMSHGF:'
        WRITE(6,*)'YW,YPW:',YW,YPW
      ENDIF

      IF (A0SCALE.NE.0.0D0) THEN

        CALL MYBFELD(XIN,YIN,ZIN,BXI,BYI,BZI,AXI,AYI,AZI)

      ELSE

        AXI=0.D0
        AYI=0.D0
        AZI=0.D0

      ENDIF

      AXR=AXI*EWS(1)
     &  +AYI*EWS(2)+AZI*EWS(3)
      AYR=AXI*EWY(1)
     &  +AYI*EWY(2)+AZI*EWY(3)
      AZR=AXI*EWZ(1)
     &  +AYI*EWZ(2)+AZI*EWZ(3)

C--- KANONISCHE VARIABLEN

      QXI= ZW
      QYI= YW
      RNENN=SQRT(1.0D0 + ZPW**2 + YPW**2)
c18.1.2018{
c      PXI = (AZR-AZR0)/BRHOABS + ZPW/RNENN
c      PYI = (AYR-AYR0)/BRHOABS + YPW/RNENN
c18.1.2018}
c18.1.2018{
      if (icharge.gt.0) then
        PXI=AZR/BRHOABS + ZPW/RNENN-pxi0
        PYI=AYR/BRHOABS + YPW/RNENN-pyi0
      else
        PXI=-AZR/BRHOABS + ZPW/RNENN-pxi0
        PYI=-AYR/BRHOABS + YPW/RNENN-pyi0
      endif
c18.1.2018}

C****************************************************************

C--- BERECHNE DIE SCHAETZWERTE MITTELS DER LINEAREN TRANFERMATRIX

C      XF =  TRANSFM(1,1)*X + TRANSFM(1,2)*XP
C     &   +  TRANSFM(1,3)*Z + TRANSFM(1,4)*ZP
C      XPF =  TRANSFM(2,1)*X + TRANSFM(2,2)*XP
C     &   +  TRANSFM(2,3)*Z + TRANSFM(2,4)*ZP
C      ZF =  TRANSFM(3,1)*X + TRANSFM(3,2)*XP
C     &   +  TRANSFM(3,3)*Z + TRANSFM(3,4)*ZP
C      ZPF =  TRANSFM(4,1)*X + TRANSFM(4,2)*XP
C     &   +  TRANSFM(4,3)*Z + TRANSFM(4,4)*ZP

C         RNENN = DSQRT(1.D0 + XPF**2 + ZPF**2)
C
C         PX =  XPF/RNENN
C         PY =  ZPF/RNENN
C

      PX =  TRANSFM(2,1)*QXI + TRANSFM(2,2)*PXI
     &  +  TRANSFM(2,3)*QYI + TRANSFM(2,4)*PYI
      PY =  TRANSFM(4,1)*QXI + TRANSFM(4,2)*PXI
     &  +  TRANSFM(4,3)*QYI + TRANSFM(4,4)*PYI

C****************************************************************

C      goto 321

c--  Beginn der Newton Fit Routine:

      DO JLOOP=1,JMAX

C****************************************************************

C---  POTENZTERME BERECHEN
C  Z.B. QXPOW(1)=0.,QXPOW(2)=1.,QXPOW(3)=QXI,QXPOW(4)=QXI**2

        DO IPOW=2,MORDNG
          QXPOW(IPOW+1)=QXPOW(IPOW)*QXI
          PXPOW(IPOW+1)=PXPOW(IPOW)*PX
          QYPOW(IPOW+1)=QYPOW(IPOW)*QYI
          PYPOW(IPOW+1)=PYPOW(IPOW)*PY
        END DO

        FFX=0.D0
        DFFX=0.D0
        FFY=0.D0
        DFFY=0.D0

C--- PARTIELLE ABLEITUNGEN DER ERZEUGENDEN-FUNKTION BERECHEN

        DO I=0,MORDNG-1
          DO J=0,MORDNG-1
            DO K=0,MORDNG-1
              DO L=0,MORDNG-1

                I1=I+1 !IN DEN FOLGENDEN BERECHNUNGEN SIND
                       !DIE I,J,K,L MATH.
                J1=J+1 !INDIZES, DIE I1,J1,K1,L1 DIE
                       ! ENTSPRECHENDEN FORTRAN
                K1=K+1 !INDIZES
                L1=L+1

                I2=I+2
                J2=J+2
                K2=K+2
                L2=L+2

                IF(I+J+K+L.GT.0 .AND. I+J+K+L.LT.MORDNG)
     &              THEN !OHNE CLOSED ORBIT

                  FFX =
     &              FFX + dble(I)*  AKOEFF(I1,J1,K1,L1)*
     &              QXPOW(I2-1)*PXPOW(J2)* QYPOW(K2)*  PYPOW(L2)

                  DFFX =
     &              DFFX + dble(I*J)*AKOEFF(I1,J1,K1,L1)*
     &              QXPOW(I2-1)*PXPOW(J2-1)*QYPOW(K2)*PYPOW(L2)
                  FFY =
     &              FFY  + dble(K)*  AKOEFF(I1,J1,K1,L1)*
     &              QXPOW(I2)*  PXPOW(J2)* QYPOW(K2-1)*  PYPOW(L2)

                  DFFY =
     &              DFFY + dble(K*L)*AKOEFF(I1,J1,K1,L1)*
     &              QXPOW(I2)*PXPOW(J2)*QYPOW(K2-1)*PYPOW(L2-1)

                ENDIF

              ENDDO
            ENDDO
          ENDDO
        ENDDO

        FFX=FFX-PXI
        FFY=FFY-PYI

C--- NEUE SCHAETZWERTE BERECHEN

        DXX=0.0
        IF(DFFX.NE.0.0) DXX  = FFX/DFFX
        DYY=0.0
        IF(DFFY.NE.0.0) DYY  = FFY/DFFY
        PX   = PX - DXX
        PY   = PY - DYY

        IF((DABS(DYY)+DABS(DXX)).LT.YACC) GOTO 123

      ENDDO

123   CONTINUE

      IF (JLOOP.GE.JMAX) THEN

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '*** WARNING SR IDTRMSHGF: NEWTOWN-FIT FAILED ***'
        WRITE(LUNGFO,*)'Z,ZP,Y,YP:'
        WRITE(LUNGFO,*)ZW,ZPW,YW,YPW
        WRITE(LUNGFO,*)'FFX,DXX,FFY,DYY:'
        WRITE(LUNGFO,*)FFX,DXX,FFY,DYY
        WRITE(LUNGFO,*)'FFX,DXX,FFY,DYY:'
        WRITE(LUNGFO,*)FFX,DXX,FFY,DYY
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)
     &    '*** WARNING SR IDTRMSHGF: NEWTOWN-FIT FAILED ***'
        WRITE(LUNGFO,*)'Z,ZP,Y,YP:'
        WRITE(LUNGFO,*)ZW,ZPW,YW,YPW
        WRITE(6,*)'FFX,DXX,FFY,DYY:'
        WRITE(6,*)FFX,DXX,FFY,DYY

      ENDIF

c-- Ende der Newton Fit Routine:

321   CONTINUE

CC Hier werden die End-Koordinaten berechnet

C****************************************************************
C---  POTENZTERME BERECHEN
C  Z.B. QXPOW(1)=0.,QXPOW(2)=1.,QXPOW(3)=QXI,QXPOW(4)=QXI**2

      DO IPOW=2,MORDNG
C          QXPOW(IPOW+1)=QXPOW(IPOW)*QXI
        PXPOW(IPOW+1)=PXPOW(IPOW)*PX
C          QYPOW(IPOW+1)=QYPOW(IPOW)*QYI
        PYPOW(IPOW+1)=PYPOW(IPOW)*PY
      END DO

C--- KOORDINATEN BERECHNEN

      QX=0.0D0
      QY=0.0D0

      DO I=0,MORDNG-1
        DO J=0,MORDNG-1
          DO K=0,MORDNG-1
            DO L=0,MORDNG-1

              I1=I+1 !IN DEN FOLGENDEN BERECHNUNGEN
                     !SIND DIE I,J,K,L MATH.
              J1=J+1 !INDIZES, DIE I1,J1,K1,L1 DIE
                     !ENTSPRECHENDEN FORTRAN
              K1=K+1 !INDIZES
              L1=L+1

              I2=I+2
              J2=J+2
              K2=K+2
              L2=L+2

              IF(I+J+K+L.GT.0 .AND. I+J+K+L.LT.MORDNG)
     &            THEN !OHNE CLOSED ORBIT

                QX =  QX + dble(J)*  AKOEFF(I1,J1,K1,L1)*
     &            QXPOW(I2)*PXPOW(J2-1)* QYPOW(K2)*  PYPOW(L2)

                QY = QY + dble(L)*  AKOEFF(I1,J1,K1,L1)*
     &            QXPOW(I2)*  PXPOW(J2)* QYPOW(K2)*  PYPOW(L2-1)

              ENDIF

            ENDDO
          ENDDO
        ENDDO
      ENDDO

C****************************************************************

1234  CONTINUE

C     UMSETZTEN DER KANONISCHEN VARIABLEN IN KOORDINATEN:

      ZW  = QX
      YW  = QY

      XL=XF0+ZW*EWZF(1)+YW*EWYF(1)
      YL=YF0+ZW*EWZF(2)+YW*EWYF(2)
      ZL=ZF0+ZW*EWZF(3)+YW*EWYF(3)

      IF (A0SCALE.NE.0.0D0) THEN

        CALL MYBFELD(XL,YL,ZL,BXF,BYF,BZF,AXF,AYF,AZF)

      ELSE

        AXF=0.D0
        AYF=0.D0
        AZF=0.D0

      ENDIF


      AXR=AXF*EWSF(1)
     &  +AYF*EWSF(2)+AZF*EWSF(3)
      AYR=AXF*EWYF(1)
     &  +AYF*EWYF(2)+AZF*EWYF(3)
      AZR=AXF*EWZF(1)
     &  +AYF*EWZF(2)+AZF*EWZF(3)

c18.1.2018{
c        PAX=PX-(AZR-AZRF0)/BRHOABS
c        PAY=PY-(AYR-AYRF0)/BRHOABS

      if (icharge.gt.0) then
        pax=px-azr/brhoabs+pxf0
        pay=py-ayr/brhoabs+pyf0
      else
        pax=px+azr/brhoabs+pxf0
        pay=py+ayr/brhoabs+pyf0
      endif
c18.1.2018}

      SNENN2 = 1.0D0 - PAX*PAX - PAY*PAY

c spezielle Einschub falls SNENN2<0 wird. Dann soll diese Subroutine
c nicht hier aussteigen sondern irgenwoanders im BETA-CODE. Dazu wird
c XPF und ZPF mit Faktor 1.D6 multipliziert:

      IF ( SNENN2 .GT. 0.0) THEN
        SNENN = DSQRT(SNENN2)
        ZPW = PAX/SNENN
        YPW = PAY/SNENN
      ELSE
        WRITE(6,*)'*** WARNING IN IDTRMSHGF:'
        WRITE(6,*)
     &    '*** NEGATIVE ROOT OF CANONICAL MOMENTUM -> INSTABILITY'
        ZPW = PX*1.0D6
        YPW = PY*1.0D6
      ENDIF

      YF=YL
      ZF=ZL

      YPF=YPW
      ZPF=ZPW

cmsh20171128 VXF=V0*EWSF(1)
cmsh20171128 VYF=VXF*YPF
cmsh20171128 VZF=VXF*ZPF

c18.1.2018{
c        VXF=V0*vxinn
c        VYF=v0*vyinn+vxf*ypf
c        VZF=v0*vzinn+vxf*zpf
      VXF=V0/sqrt(1.0d0+zpf**2+ypf**2)
      VYF=VXF*YPF
      VZF=VXF*ZPF
c18.1.2018}

      RETURN
      END
