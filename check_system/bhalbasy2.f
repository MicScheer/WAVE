*CMZ :  4.00/15 24/03/2022  07.14.23  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.15.16  by  Michael Scheer
*CMZ :  3.05/04 28/06/2018  08.21.29  by  Michael Scheer
*CMZ :  3.04/00 16/01/2018  10.03.31  by  Michael Scheer
*CMZ :  3.03/02 17/12/2015  10.12.49  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  11.18.26  by  Michael Scheer
*CMZ :  2.62/04 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.20/09 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.03.48  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.45.16  by  Michael Scheer
*CMZ : 00.01/12 27/09/96  15.28.21  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.23  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.59  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BHALBASY2(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)

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

C SUBROUTINE CALCULATES MAGNETIC FIELD AND VECTOR POTENTIAL FOR SIMPLE
C WAVELENGTH SHIFTER MODEL WITH END POLES. THE FIELD OF THE SINGLE POLES
C CORRESPONDS TO HALBACHS FORMULA.
C INPUT AND OUTPUT CORRESPOND TO LAB.-SYSTEM, WHERE X IS COORDINATE ON
C LONGITUDINAL AXIS
C IF FLAG IAHWFOUR IS SET FIELD IS SUPERPOSITION OF HALBACH-WIGGLER
C ACCORDING TO FOURIER EXPANSION OF ON-AXIS FIELD

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,halbasy.
      include 'halbasy.cmn'
*KEND.

      INTEGER ICAL,NAHWPOL

      DOUBLE PRECISION PARK,WLEN1,EHARM1

      DOUBLE PRECISION XKHALBASY2,YKHALBASY2,ZKHALBASY2,ZLHALBASY2,X2

      DOUBLE PRECISION XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT,
     &         XKX,YKY,ZKZ,DSNXKX,DCSXKX,DSHYKY,DCHYKY,DSNZKZ,DCSZKZ
     &        ,BXH,BYH,BZH,AXH,AYH,AZH,AHWMOD,X,TOTLEN,TOTLEN2

*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DATA TOTLEN/1.D30/
      DATA TOTLEN2/1.D30/

      DATA ICAL/0/

C--- K-VALUES

      XKHALBASY=0.0D0
      YKHALBASY=0.0D0
      ZKHALBASY=0.0D0

      IF (ZLHALBASY.NE.0.0D0) ZKHALBASY=2.D0*PI1/ZLHALBASY
      IF (YLHALBASY.NE.0.0D0) YKHALBASY=2.D0*PI1/YLHALBASY
      IF (XLHALBASY.NE.0.0D0) XKHALBASY=2.D0*PI1/XLHALBASY

      IF (IAHWFOUR.NE.0 .AND. XLHALBASY.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN BHALBASY2 ***'
        WRITE(LUNGFO,*)
     &    'IF IAHWFOUR IS SET, XLHALBASY MUST BE 0.0'
        WRITE(LUNGFO,*)
     &    'SPECIFY THE GRADIENT VIA XLENFOUR IN NAMELIST FOURIER'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN BHALBASY2 ***'
        WRITE(6,*)
     &    'IF IAHWFOUR IS SET, XLHALBASY MUST BE 0.0'
        WRITE(6,*)
     &    'THE GRADIENT IS GIVEN BY XLENFOUR IN NAMELIST FOURIER'
        WRITE(6,*)
        STOP
      ENDIF

C--- ADJUST K-VALUES

      YKHALBASY=DSQRT(ZKHALBASY**2+XKHALBASY**2)
      YLHALBASY=2.D0*PI1/YKHALBASY

      IF (KHALBASY.NE.0.OR.IAHWFOUR.NE.0) THEN

C--- BENDING RADIUS AND DEVICE LENGTH

        IF(B0HALBASY.NE.0.00D0) THEN
          RHALBASY=DMYGAMMA*EMASSE1/(CLIGHT1*B0HALBASY)
        ELSE
          RHALBASY=0.00D0
        ENDIF

        PARK=ECHARGE1*DABS(B0HALBASY)*ZLHALBASY/(2.*PI1*EMASSKG1*CLIGHT1)
        WLEN1=(1+PARK**2/2.)/2./DMYGAMMA**2*ZLHALBASY*1.D9

        TOTLEN=ZLHALBASY*((AHWPOL-1.D0)/2.D0+1.D0)
        TOTLEN2=TOTLEN/2.D0

        IF (WLEN1.NE.0.00) EHARM1=WTOE1/WLEN1

        IF (ICAL.EQ.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '     Parameters of simple wavelength shifter model:'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '     peak field [T] and bending radius [m]:  ',
     &      SNGL(B0HALBASY),SNGL(RHALBASY)
          IF (IAHWFOUR.EQ.0) THEN
            WRITE(LUNGFO,*)
     &        '     l0, l0x, l0y [m]: ',
     &        SNGL(ZLHALBASY),SNGL(XLHALBASY),SNGL(YLHALBASY)
            WRITE(LUNGFO,*)
     &        '     k, kx, ky [1/m]:  ',
     &        SNGL(ZKHALBASY),SNGL(XKHALBASY),SNGL(YKHALBASY)
            WRITE(LUNGFO,*)
          ELSE
            WRITE(LUNGFO,*)
     &        '     l0, l0y [m]: ',
     &        SNGL(ZLHALBASY),SNGL(YLHALBASY)
            WRITE(LUNGFO,*)
     &        '     k, ky [1/m]:  ',
     &        SNGL(ZKHALBASY),SNGL(YKHALBASY)
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
     &        '     *** flag IAHWFOUR is set, i.e. global kx is taken'
            WRITE(LUNGFO,*)
     &        '     from XLENFOUR of namelist FOURIER.'
          ENDIF
          WRITE(LUNGFO,*)
     &      '     peak field ratio, total number of poles:',
     &      SNGL(FASYM),NINT(AHWPOL+2.D0)
          WRITE(LUNGFO,*)
     &      '     total device length, half length:       ',
     &      SNGL(TOTLEN),SNGL(TOTLEN2)
          WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)
     &       '     x-position of device center:            ',
     &       SNGL(xcenhal)
          WRITE(LUNGFO,*)
     &      '      Deflection parameter K, 1. harmonical [eV] (main poles):'
          WRITE(LUNGFO,*)
     &      '     ',SNGL(PARK),SNGL(EHARM1)
          WRITE(LUNGFO,*)'      Critical energy [eV] (main poles):'
     &      ,SNGL(ecdipev1*DABS(B0HALBASY)*DMYENERGY**2)
          WRITE(LUNGFO,*)
        ENDIF

        ICAL=1
      ENDIF

      TOTLEN=ZLHALBASY*((AHWPOL-1.D0)/2.D0+1.D0)
      TOTLEN2=TOTLEN/2.D0

      IF(IAHWFOUR.NE.0) THEN
          CALL BFOUR(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)
          RETURN
      ENDIF

      NAHWPOL=AHWPOL
      AHWMOD=-ISIGN(1,-(MOD(NAHWPOL,4)-2))/2.D0

      X=XIN
      IF (DABS(XIN).GT.TOTLEN2) THEN
          BXOUT=0.0
          BYOUT=0.0
          BZOUT=0.0
          AXOUT=0.0
          AYOUT=0.0
          AZOUT=0.0
          RETURN
      ENDIF

      IF (DABS(X).LE.TOTLEN2-ZLHALBASY/2.D0) THEN

          XKX=XKHALBASY*(-ZIN)
          YKY=YKHALBASY*YIN
          ZKZ=ZKHALBASY*X

          DSNXKX=DSIN(XKX)
          DCSXKX=DCOS(XKX)
          DSHYKY=DSINH(YKY)
          DCHYKY=DSQRT(1.D0+DSHYKY*DSHYKY)
          DSNZKZ=DSIN(ZKZ)
          DCSZKZ=DCOS(ZKZ)


          BXH=-XKHALBASY/YKHALBASY*B0HALBASY*DSNXKX*DSHYKY*DCSZKZ
          BYH=                 B0HALBASY*DCSXKX*DCHYKY*DCSZKZ
          BZH=-ZKHALBASY/YKHALBASY*B0HALBASY*DCSXKX*DSHYKY*DSNZKZ

          AXH=B0HALBASY/ZKHALBASY*                    DCSXKX*DCHYKY*DSNZKZ
          AYH=B0HALBASY/ZKHALBASY*XKHALBASY/YKHALBASY*DSNXKX*DSHYKY*DSNZKZ
          AZH=0.0

          BZOUT=-BXH
          BYOUT=BYH
          BXOUT=BZH

          AZOUT=-AXH
          AYOUT=AYH
          AXOUT=AZH

          RETURN

      ELSE

        XKHALBASY2=XKHALBASY

        ZKHALBASY2=ZKHALBASY
        ZLHALBASY2=2.D0*PI1/ZKHALBASY2
          YKHALBASY2=DSQRT(ZKHALBASY2**2+XKHALBASY2**2)
          YLHALBASY=2.D0*PI1/YKHALBASY2

        X2=X+TOTLEN2+ZLHALBASY/2.D0

          XKX=XKHALBASY2*(-ZIN)
          YKY=YKHALBASY2*YIN
          ZKZ=ZKHALBASY2*(X2)

          DSNXKX=DSIN(XKX)
          DCSXKX=DCOS(XKX)
          DSHYKY=DSINH(YKY)
          DCHYKY=DSQRT(1.D0+DSHYKY*DSHYKY)
          DSNZKZ=DSIN(ZKZ)
          DCSZKZ=DCOS(ZKZ)

          BXH=-XKHALBASY2/YKHALBASY2*B0HALBASY*DSNXKX*DSHYKY*DCSZKZ
          BYH=                 B0HALBASY*DCSXKX*DCHYKY*DCSZKZ
          BZH=-ZKHALBASY2/YKHALBASY2*B0HALBASY*DCSXKX*DSHYKY*DSNZKZ

          AXH=B0HALBASY/ZKHALBASY2*                    DCSXKX*DCHYKY*DSNZKZ
          AYH=B0HALBASY/ZKHALBASY2*XKHALBASY2/YKHALBASY2*DSNXKX*DSHYKY*DSNZKZ
          AZH=0.0

          ZKHALBASY2=ZKHALBASY*2.D0
          ZLHALBASY2=2.D0*PI1/ZKHALBASY2
          YKHALBASY2=DSQRT(ZKHALBASY2**2+XKHALBASY2**2)
          YLHALBASY=2.D0*PI1/YKHALBASY2

          XKX=XKHALBASY2*(-ZIN)
          YKY=YKHALBASY2*YIN
          ZKZ=ZKHALBASY2*(X2)

          DSNXKX=DSIN(XKX)
          DCSXKX=DCOS(XKX)
          DSHYKY=DSINH(YKY)
          DCHYKY=DSQRT(1.D0+DSHYKY*DSHYKY)
          DSNZKZ=DSIN(ZKZ)
          DCSZKZ=DCOS(ZKZ)

          BXH=BXH-XKHALBASY2/YKHALBASY2*B0HALBASY*DSNXKX*DSHYKY*DCSZKZ
          BYH=BYH+                      B0HALBASY*DCSXKX*DCHYKY*DCSZKZ
          BZH=BZH-ZKHALBASY2/YKHALBASY2*B0HALBASY*DCSXKX*DSHYKY*DSNZKZ

          AXH=AXH+B0HALBASY/ZKHALBASY2*                    DCSXKX*DCHYKY*DSNZKZ
          AYH=AYH+B0HALBASY/ZKHALBASY2*XKHALBASY2/YKHALBASY2*DSNXKX*DSHYKY*DSNZKZ
          AZH=0.0

          BZOUT=BXH*AHWMOD
          BYOUT=-BYH*AHWMOD
          BXOUT=-BZH*AHWMOD

          AZOUT=AXH*AHWMOD
          AYOUT=-AYH*AHWMOD
          AXOUT=-AZH*AHWMOD

      RETURN

      ENDIF

      END
