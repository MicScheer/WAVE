*CMZ :  4.00/11 07/05/2021  07.29.16  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.15.16  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  11.18.26  by  Michael Scheer
*CMZ :  2.64/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.36/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ : 00.01/10 21/08/96  12.30.33  by  Michael Scheer
*CMZ : 00.01/04 30/11/94  14.09.41  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.18  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.57  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BHALBAstandalone(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT
     &  ,AXOUT,AYOUT,AZOUT,b0halba,xlhalba,ylhalba,zlhalba,perhal)

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

C Calculates magnetic field according to Halbach's formulas
C Inside the routine the coodinate system corresponds to
C Halbach's convention i.e. Z is the longitudinal coordinate
C Input and Output are converted according to the coordinate
C system of WAVE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEND.

*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER ICAL

      DOUBLE PRECISION XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT,
     &  XKX,YKY,ZKZ,DSNXKX,DCSXKX,DSHYKY,DCHYKY,DSNZKZ,DCSZKZ
     &  ,BXH,BYH,BZH,AXH,AYH,AZH,perhal,b0halba,xkhalba,ykhalba,zkhalba,
     &  zlenhal2,zlhalba,xlhalba,ylhalba,zlenhal

      DOUBLE PRECISION PARK,WLEN1,EHARM1

      save

      DATA ICAL/0/

      IF (ICAL.EQ.0) THEN

        IF (PERHAL.LT.0.0d0) THEN
          ZLHALBA=-ZLENHAL/PERHAL
          XLHALBA= XLHALBA/PERHAL
        ENDIF

        XKHALBA=0.0d0
        YKHALBA=0.0d0
        ZKHALBA=0.0d0

        IF (ZLHALBA.NE.0.0d0) ZKHALBA=2.0d0*PI1/ZLHALBA
        IF (YLHALBA.NE.0.0d0) YKHALBA=2.0d0*PI1/YLHALBA
        IF (XLHALBA.NE.0.0d0) XKHALBA=2.0d0*PI1/XLHALBA

C--- ADJUST K-VALUES

        YKHALBA=DSQRT(ZKHALBA*ZKHALBA+XKHALBA*XKHALBA)
        YLHALBA=2.0d0*PI1/YKHALBA

        IF(PERHAL.GE.0) ZLENHAL=PERHAL*ZLHALBA     !29.10.91

        ZLENHAL2=ZLENHAL/2.0d0

C--- BENDING RADIUS

        PARK=ECHARGE1*DABS(B0HALBA)*ZLHALBA/(2.*PI1*EMASSKG1*CLIGHT1)
        WLEN1=(1+PARK**2/2.)/2./DMYGAMMA**2*ZLHALBA*1.D9
        IF (WLEN1.NE.0.0) EHARM1=WTOE1/WLEN1

        ICAL =1

      ENDIF

      IF (KHALBA.LT.0.AND.ABS(XIN).GT.ZLENHAL2) THEN
        BXOUT=0.0D0
        BYOUT=0.0D0
        BZOUT=0.0D0
        AXOUT=0.0D0
        AYOUT=0.0D0
        AZOUT=0.0D0
        RETURN
      ENDIF

      XKX=XKHALBA*(-ZIN)
      YKY=YKHALBA*YIN
      ZKZ=ZKHALBA*XIN

      DSNXKX=DSIN(XKX)
      DCSXKX=DCOS(XKX)
      DSHYKY=DSINH(YKY)
      DCHYKY=DSQRT(1.D0+DSHYKY*DSHYKY)
      DSNZKZ=DSIN(ZKZ)
      DCSZKZ=DCOS(ZKZ)

      BXH=-XKHALBA/YKHALBA*B0HALBA*DSNXKX*DSHYKY*DCSZKZ
      BYH=                 B0HALBA*DCSXKX*DCHYKY*DCSZKZ
      BZH=-ZKHALBA/YKHALBA*B0HALBA*DCSXKX*DSHYKY*DSNZKZ

      AXH=B0HALBA/ZKHALBA*                DCSXKX*DCHYKY*DSNZKZ
      AYH=B0HALBA/ZKHALBA*XKHALBA/YKHALBA*DSNXKX*DSHYKY*DSNZKZ
      AZH=0.

      BZOUT=-BXH
      BYOUT=BYH
      BXOUT=BZH

      AZOUT=-AXH
      AYOUT= AYH
      AXOUT= AZH

      RETURN
      END
