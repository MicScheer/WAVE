*CMZ :  4.00/15 19/03/2022  09.33.42  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.15.16  by  Michael Scheer
*CMZ :  3.05/04 28/06/2018  08.22.42  by  Michael Scheer
*CMZ :  3.03/02 17/12/2015  10.13.22  by  Michael Scheer
*CMZ :  3.01/00 10/04/2013  09.49.08  by  Michael Scheer
*CMZ :  3.00/02 10/04/2013  09.24.32  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  11.18.26  by  Michael Scheer
*CMZ :  2.67/02 03/05/2012  09.32.01  by  Michael Scheer
*CMZ :  2.62/04 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.52/06 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.20/09 22/03/2001  13.10.36  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.03.48  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.45.16  by  Michael Scheer
*CMZ : 00.01/12 27/09/96  15.28.21  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.23  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.59  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BHALBASY(XINI,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)

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
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,halbasy.
      include 'halbasy.cmn'
*KEND.

      INTEGER ICAL,NAHWPOL

      DOUBLE PRECISION PARK,WLEN1,EHARM1

      DOUBLE PRECISION XIN,XINI,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT,
     &         XKX,YKY,ZKZ,DSNXKX,DCSXKX,DSHYKY,DCHYKY,DSNZKZ,DCSZKZ
     &        ,BXH,BYH,BZH,AXH,AYH,AZH,AHWMOD,X,TOTLEN,TOTLEN2

*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DATA TOTLEN/1.0d30/
      DATA TOTLEN2/1.0d30/

      DATA ICAL/0/

      IF (ICAL.EQ.0) THEN

        park=pkhalbasy

        if (nhhalbasy.ne.0.and.hhalbasy.ne.0.0d0) then
          if (hhalbasy.eq.-9999.0d0) then
            if (ifreq2p.eq.1) then
              hhalbasy=freqlow
            else
              hhalbasy=(freqlow+freqhig)/2.0d0
            endif
          else if (hhalbasy.lt.0.0d0) then
            hhalbasy=-wtoe1/hhalbasy
          endif
          WLEN1=wtoe1/abs(hhalbasy/nhhalbasy)
          park=2.0d0*(wlen1/(zlhalbasy*1.0D9/2.0d0/DMYGAMMA**2)-1.0d0)
          if (park.lt.0.0d0) then
            write(6,*)
     &        '*** Error in BHALBASY:'
            write(6,*)
     &        'Inconsistent values of NHHALBASY, HHALBASY, and ZLHALBASY'
            write(6,*)' '
            write(lungfo,*)
     &        '*** Error in BHALBASY:'
            write(lungfo,*)
     &        'Inconsistent values of NHHALBASY, HHALBASY, and ZLHALBASY'
            write(lungfo,*)' '
            stop
          endif
          park=sqrt(park)
          pkhalbasy=park
        endif

        IF (park.ne.0.0d0) THEN
          B0halbasy=pkhalbasy/(echarge1*zlhalbasy/(2.*PI1*EMASSKG1*CLIGHT1))
        endif

      endif !ical

      XIN=XINI-XCENHAL

      IF (FASYM.EQ.2.0d0)  THEN
        CALL BHALBASY2(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)
        RETURN
      ENDIF


      IF (ICAL.EQ.0) THEN

C--- K-VALUES

         XKHALBASY=0.0d0
         YKHALBASY=0.0d0
         ZKHALBASY=0.0d0

         IF (ZLHALBASY.NE.0.0d0) ZKHALBASY=2.0d0*PI1/ZLHALBASY
         IF (YLHALBASY.NE.0.0d0) YKHALBASY=2.0d0*PI1/YLHALBASY
         IF (XLHALBASY.NE.0.0d0) XKHALBASY=2.0d0*PI1/XLHALBASY

         IF (IAHWFOUR.NE.0 .AND. XLHALBASY.NE.0) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN BHALBASY ***'
            WRITE(LUNGFO,*)
     & 'IF IAHWFOUR IS SET, XLHALBASY MUST BE 0.'
            WRITE(LUNGFO,*)
     & 'SPECIFY THE GRADIENT VIA XLENFOUR IN NAMELIST FOURIER'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN BHALBASY ***'
            WRITE(6,*)
     & 'IF IAHWFOUR IS SET, XLHALBASY MUST BE 0.'
            WRITE(6,*)
     & 'THE GRADIENT IS GIVEN BY XLENFOUR IN NAMELIST FOURIER'
            WRITE(6,*)
            STOP
         ENDIF

C--- ADJUST K-VALUES

         YKHALBASY=DSQRT(ZKHALBASY**2+XKHALBASY**2)
         YLHALBASY=2.0d0*PI1/YKHALBASY

         IF (KHALBASY.NE.0.OR.IAHWFOUR.NE.0) THEN

C--- BENDING RADIUS AND DEVICE LENGTH

           IF(B0HALBASY.NE.0.0D0) THEN
             RHALBASY=DMYGAMMA*EMASSE1/(CLIGHT1*B0HALBASY)
           ELSE
             RHALBASY=0.0D0
           ENDIF

           TOTLEN=ZLHALBASY*(AHWPOL+FASYM)/2.0d0
           TOTLEN2=TOTLEN/2.0d0

           PARK=ECHARGE1*DABS(B0HALBASY)*ZLHALBASY/(2.*PI1*EMASSKG1*CLIGHT1)
           WLEN1=(1+PARK**2/2.)/2./DMYGAMMA**2*ZLHALBASY*1.0d9
           IF (WLEN1.NE.0.0) EHARM1=WTOE1/WLEN1

           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)
     &       '     Parameters of simple wavelength shifter model:'
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)
     &       '     peak field [T] and bending radius [m]:  ',
     &       SNGL(B0HALBASY),SNGL(RHALBASY)
           IF (IAHWFOUR.EQ.0) THEN
             WRITE(LUNGFO,*)
     &         '     l0, l0x, l0y [m]: ',
     &         SNGL(ZLHALBASY),SNGL(XLHALBASY),SNGL(YLHALBASY)
             WRITE(LUNGFO,*)
     &         '     k, kx, ky [1/m]:  ',
     &         SNGL(ZKHALBASY),SNGL(XKHALBASY),SNGL(YKHALBASY)
             WRITE(LUNGFO,*)
           ELSE
             WRITE(LUNGFO,*)
     &         '     l0, l0y [m]: ',
     &         SNGL(ZLHALBASY),SNGL(YLHALBASY)
             WRITE(LUNGFO,*)
     &         '     k, ky [1/m]:  ',
     &         SNGL(ZKHALBASY),SNGL(YKHALBASY)
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
     &         '     *** flag IAHWFOUR is set, i.e. global kx is taken'
             WRITE(LUNGFO,*)
     &         '     from XLENFOUR of namelist FOURIER.'
           ENDIF
           WRITE(LUNGFO,*)
     &       '     peak field ratio, total number of poles:',
     &       SNGL(FASYM),NINT(AHWPOL+2.0d0)
           WRITE(LUNGFO,*)
     &       '     total device length, half length:       ',
     &       SNGL(TOTLEN),SNGL(TOTLEN2)
           WRITE(LUNGFO,*)
     &       '     x-position of device center:            ',
     &       SNGL(xcenhal)
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)
     &       '      Deflection parameter K, 1. harmonical [eV] (main poles):'
           WRITE(LUNGFO,*)
     &       '     ',SNGL(PARK),SNGL(EHARM1)
           WRITE(LUNGFO,*)'      Critical energy [eV] (main poles):'
     &       ,SNGL(ecdipev1*DABS(B0HALBASY)*DMYENERGY**2)
           WRITE(LUNGFO,*)
         ENDIF

         ICAL=1
       ENDIF

       IF(IAHWFOUR.NE.0) THEN
         CALL BFOUR(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)
         RETURN
       ENDIF

       NAHWPOL=AHWPOL
       AHWMOD=ISIGN(1,-(MOD(NAHWPOL,4)-2))

       X=XIN
       IF (DABS(XIN).GT.TOTLEN2) THEN
         BXOUT=0.0d0
         BYOUT=0.0d0
         BZOUT=0.0d0
         AXOUT=0.0d0
         AYOUT=0.0d0
         AZOUT=0.0d0
         RETURN
       ENDIF

       IF (DABS(X).LE.(AHWPOL*ZLHALBASY)/4.) THEN

         XKX=XKHALBASY*(-ZIN)
         YKY=YKHALBASY*YIN
         ZKZ=ZKHALBASY*X

         DSNXKX=DSIN(XKX)
         DCSXKX=DCOS(XKX)
         DSHYKY=DSINH(YKY)
         DCHYKY=DSQRT(1.0d0+DSHYKY*DSHYKY)
         DSNZKZ=DSIN(ZKZ)
         DCSZKZ=DCOS(ZKZ)

         BXH=-XKHALBASY/YKHALBASY*B0HALBASY*DSNXKX*DSHYKY*DCSZKZ
         BYH=                 B0HALBASY*DCSXKX*DCHYKY*DCSZKZ
         BZH=-ZKHALBASY/YKHALBASY*B0HALBASY*DCSXKX*DSHYKY*DSNZKZ

         AXH=B0HALBASY/ZKHALBASY*                    DCSXKX*DCHYKY*DSNZKZ
         AYH=B0HALBASY/ZKHALBASY*XKHALBASY/YKHALBASY*DSNXKX*DSHYKY*DSNZKZ
         AZH=0.0d0

         BZOUT=-BXH
         BYOUT=BYH
         BXOUT=BZH

         AZOUT=-AXH
         AYOUT=AYH
         AXOUT=AZH

         RETURN

       ELSE IF (X.GT.(AHWPOL*ZLHALBASY)/4.) THEN

         XKX=XKHALBASY*(-ZIN/FASYM*2.)
         YKY=YKHALBASY*YIN/FASYM*2.
         ZKZ=(X*ZKHALBASY-AHWPOL*PI1/2.)/FASYM*2.-PI1/2.

         DSNXKX=DSIN(XKX)
         DCSXKX=DCOS(XKX)
         DSHYKY=DSINH(YKY)
         DCHYKY=DSQRT(1.0d0+DSHYKY*DSHYKY)
         DSNZKZ=DSIN(ZKZ)
         DCSZKZ=DCOS(ZKZ)

         BXH=-XKHALBASY/YKHALBASY*B0HALBASY/(-FASYM)*DSNXKX*DSHYKY*DCSZKZ
         BYH=                 B0HALBASY/(-FASYM)*DCSXKX*DCHYKY*DCSZKZ
         BZH=-ZKHALBASY/YKHALBASY*B0HALBASY/(-FASYM)*DCSXKX*DSHYKY*DSNZKZ

         AXH=B0HALBASY/(-FASYM)/ZKHALBASY*           DCSXKX*DCHYKY*DSNZKZ
         AYH=B0HALBASY/(-FASYM)/ZKHALBASY*XKHALBASY/YKHALBASY
     &     *DSNXKX*DSHYKY*DSNZKZ
         AZH=0.0d0

         BZOUT=-BXH*AHWMOD
         BYOUT=BYH*AHWMOD
         BXOUT=BZH*AHWMOD

         AZOUT=-AXH*AHWMOD
         AYOUT=AYH*AHWMOD
         AXOUT=AZH*AHWMOD

         RETURN

       ELSE

         XKX=XKHALBASY*(-ZIN/FASYM*2.)
         YKY=YKHALBASY*YIN/FASYM*2.
         ZKZ=((X*ZKHALBASY+AHWPOL*PI1/2.)/FASYM*2.+PI1/2.)

         DSNXKX=DSIN(XKX)
         DCSXKX=DCOS(XKX)
         DSHYKY=DSINH(YKY)
         DCHYKY=DSQRT(1.0d0+DSHYKY*DSHYKY)
         DSNZKZ=DSIN(ZKZ)
         DCSZKZ=DCOS(ZKZ)

         BXH=-XKHALBASY/YKHALBASY*B0HALBASY/(-FASYM)*DSNXKX*DSHYKY*DCSZKZ
         BYH=                 B0HALBASY/(-FASYM)*DCSXKX*DCHYKY*DCSZKZ
         BZH=-ZKHALBASY/YKHALBASY*B0HALBASY/(-FASYM)*DCSXKX*DSHYKY*DSNZKZ

         AXH=B0HALBASY/(-FASYM)/ZKHALBASY*           DCSXKX*DCHYKY*DSNZKZ
         AYH=B0HALBASY/(-FASYM)/ZKHALBASY*XKHALBASY/YKHALBASY
     &     *DSNXKX*DSHYKY*DSNZKZ
         AZH=0.0d0

         BZOUT=-BXH*AHWMOD
         BYOUT=BYH*AHWMOD
         BXOUT=BZH*AHWMOD

         AZOUT=-AXH*AHWMOD
         AYOUT=AYH*AHWMOD
         AXOUT=AZH*AHWMOD

         RETURN

       ENDIF

      END
