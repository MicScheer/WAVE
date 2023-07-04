*CMZ :  3.00/01 26/03/2013  19.46.50  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.70/06 04/01/2013  13.25.02  by  Michael Scheer
*CMZ :  2.67/02 19/04/2012  09.50.56  by  Michael Scheer
*CMZ :  2.20/09 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.13/08 02/03/2000  11.25.34  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.35  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.43.01  by  Michael Scheer
*CMZ : 00.01/06 15/02/95  20.20.16  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.32.24  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.48.56  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.39  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BUCROSS_AX(XINOLD,YIN,ZIN,BXOUT,BYOUT,BZOUT,
     &                               AXOUT,AYOUT,AZOUT)
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


C--- MAGNETIC FIELD OF A CROSSED UNDULATOR. FIRST UNDULATOR IS HORIZONTAL,
C    MODULATOR IS HORIZONTAL, SECOND UNDULATOR IS VERTICAL
C    MODULATOR CENTER IS DEVICE CENTER

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,halbasy.
      include 'halbasy.cmn'
*KEEP,ucross.
      include 'ucross.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER ICAL,I

      DOUBLE PRECISION X,Y,Z,BY,BZ,BXOUT,BYOUT,BZOUT,AY,AZ,AXOUT,AYOUT,AZOUT
      DOUBLE PRECISION XIN,YIN,ZIN

      DOUBLE PRECISION WLEN1,EHARM1(3),PARK(3)
      DOUBLE PRECISION XCEN,XINOLD

      double precision b0halbasyo,xkhalbasyo,zlhalbasyo,zkhalbasyo,fasymo,
     &  ahwpolo,xcenhalo

      DATA ICAL/0/

      b0halbasyo=b0halbasy
      xkhalbasyo=xkhalbasy
      zlhalbasyo=zlhalbasy
      zkhalbasyo=zkhalbasy
      fasymo=fasym
      ahwpolo=ahwpol
      xcenhalo=xcenhal

      XIN=XINOLD

C--- FIRST CALL

      IF (ICAL.EQ.0) THEN


          DO I=1,3

            IF (UASYM(I).EQ.2.D0) UASYM(I)=2.00000001D0

               ZKUHAL(I)=2.*PI1/ZLUHAL(I)

              IF (XLUHAL(I).NE.0.D0) THEN
                XKUHAL(I)=2.*PI1/XLUHAL(I)
              ELSE
                XKUHAL(I)=0.0D0
              ENDIF

            YKUHAL(I)=SQRT(XKUHAL(I)*XKUHAL(I)+ZKUHAL(I)*ZKUHAL(I))
            YLUHAL(I)=2.*PI1/YKUHAL(I)

              PARK(I)=
     &          ECHARGE1*DABS(B0UCROSS(I))*ZLUHAL(I)
     &          /(2.*PI1*EMASSKG1*CLIGHT1)
              WLEN1=(1+PARK(I)**2/2.)/2./DMYGAMMA**2*ZLUHAL(I)*1.D9
              IF (WLEN1.NE.0.0) EHARM1(I)=WTOE1/WLEN1
              ucharm1=eharm1

          ENDDO   !I

          ULIMI(2)=-ZLUHAL(2)/2.*(DFLOAT(NMUPOL(2))+UASYM(2))/2.
          ULIMI(3)=-ULIMI(2)
          ULIMI(1)=ULIMI(2)-ZLUHAL(1)/2.*(DFLOAT(NMUPOL(1))+UASYM(1))
          ULIMI(4)=ULIMI(3)+ZLUHAL(3)/2.*(DFLOAT(NMUPOL(3))+UASYM(3))

          IF (IUCRSAX.EQ.2) THEN
             XCEN=(ULIMI(2)+ULIMI(1))/2.D0
          ELSE IF (IUCRSAX.EQ.3) THEN
             XCEN=(ULIMI(3)+ULIMI(4))/2.D0
          ELSE
             XCEN=0.0D0
          ENDIF

          ULIMI(1)=ULIMI(1)-XCEN
          ULIMI(2)=ULIMI(2)-XCEN
          ULIMI(3)=ULIMI(3)-XCEN
          ULIMI(4)=ULIMI(4)-XCEN

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '     SR BUCROSS_AX: Parameter of crossed undulator:'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     Peak field B0:    ',(SNGL(B0UCROSS(I)),I=1,3)
          WRITE(LUNGFO,*)'     Number of main poles:',(NMUPOL(I),I=1,3)
          WRITE(LUNGFO,*)'     asymmetrie of end poles:',(UASYM(I),I=1,3)
          WRITE(LUNGFO,*)'     Period lengths:   ',(SNGL(ZLUHAL(I)),I=1,3)
          WRITE(LUNGFO,*)'     2*PI/Kx:          ',(SNGL(XLUHAL(I)),I=1,3)
          WRITE(LUNGFO,*)'     Limits of modules:'
          WRITE(LUNGFO,*)'     ',(SNGL(ULIMI(I)),I=1,4)
          WRITE(LUNGFO,*)'     Deflection parameter K:',(SNGL(PARK(I)),I=1,3)
          WRITE(LUNGFO,*)
     &      '     First harmonical [eV]: ',(SNGL(EHARM1(I)),I=1,3)
          WRITE(LUNGFO,*)

          ICAL=1

      ENDIF

      IF (XIN.LT.ULIMI(1).OR.XIN.GT.ULIMI(4))THEN

          BXOUT=0.
          BYOUT=0.
          BZOUT=0.
          AXOUT=0.
          AYOUT=0.
          AZOUT=0.

      ELSE IF(XIN.LT.ULIMI(2)) THEN

C HORIZONTAL

      B0HALBASY=B0UCROSS(1)
      XKHALBASY=XKUHAL(1)
      XLHALBASY=XLUHAL(1)
      ZLHALBASY=ZLUHAL(1)
      ZKHALBASY=ZKUHAL(1)
      YKHALBASY=YKUHAL(1)
      YLHALBASY=YLUHAL(1)
      FASYM=UASYM(1)
      AHWPOL=DFLOAT(NMUPOL(1))

      X=XIN-ULIMI(1)-(ULIMI(2)-ULIMI(1))/2.
      CALL BHALBASY(X,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)

      IF     (XIN.LT.ULIMI(1)+ZLUHAL(1)*FASYM/4.) THEN
          BXOUT=BXOUT*0.5
          BYOUT=BYOUT*0.5
          BZOUT=BZOUT*0.5
      ELSE IF (XIN.LT.ULIMI(1)+ZLUHAL(1)*(FASYM/4.+1./2.)) THEN
          BXOUT=BXOUT*0.75
          BYOUT=BYOUT*0.75
          BZOUT=BZOUT*0.75
      ELSE IF (XIN.GT.ULIMI(2)-ZLUHAL(1)*FASYM/4.) THEN
          BXOUT=BXOUT*0.5
          BYOUT=BYOUT*0.5
          BZOUT=BZOUT*0.5
      ELSE IF (XIN.GT.ULIMI(2)-ZLUHAL(1)*(FASYM/4.+1./2.)) THEN
          BXOUT=BXOUT*0.75
          BYOUT=BYOUT*0.75
          BZOUT=BZOUT*0.75
      ENDIF

      ELSE IF(XIN.LT.ULIMI(3)) THEN

C MODULATOR

      B0HALBASY=B0UCROSS(2)
      XKHALBASY=XKUHAL(2)
      XLHALBASY=XLUHAL(2)
      ZLHALBASY=ZLUHAL(2)
      ZKHALBASY=ZKUHAL(2)
      YKHALBASY=YKUHAL(2)
      YLHALBASY=YLUHAL(2)
      FASYM=UASYM(2)
      AHWPOL=DFLOAT(NMUPOL(2))

      X=XIN+XCEN
      CALL BHALBASY(X,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)

      ELSE IF(XIN.LT.ULIMI(4)) THEN

C VERTICAL


      B0HALBASY=B0UCROSS(3)
      XKHALBASY=XKUHAL(3)
      XLHALBASY=XLUHAL(3)
      ZLHALBASY=ZLUHAL(3)
      ZKHALBASY=ZKUHAL(3)
      YKHALBASY=YKUHAL(3)
      YLHALBASY=YLUHAL(3)
      FASYM=UASYM(3)
      AHWPOL=DFLOAT(NMUPOL(3))

      X=XIN-ULIMI(3)-(ULIMI(4)-ULIMI(3))/2.
      Z=YIN
      Y=ZIN

      CALL BHALBASY(X,Y,Z,BXOUT,BY,BZ,AXOUT,AY,AZ)

      BZOUT= BY
      BYOUT= BZ
      AZOUT= AY
      AYOUT= AZ

      IF     (XIN.LT.ULIMI(3)+ZLUHAL(3)*FASYM/4.) THEN
          BXOUT=BXOUT*0.5
          BYOUT=BYOUT*0.5
          BZOUT=BZOUT*0.5
      ELSE IF (XIN.LT.ULIMI(3)+ZLUHAL(3)*(FASYM/4.+1./2.)) THEN
          BXOUT=BXOUT*0.75
          BYOUT=BYOUT*0.75
          BZOUT=BZOUT*0.75
      ELSE IF (XIN.GT.ULIMI(4)-ZLUHAL(3)*FASYM/4.) THEN
          BXOUT=BXOUT*0.5
          BYOUT=BYOUT*0.5
          BZOUT=BZOUT*0.5
      ELSE IF (XIN.GT.ULIMI(4)-ZLUHAL(3)*(FASYM/4.+1./2.)) THEN
          BXOUT=BXOUT*0.75
          BYOUT=BYOUT*0.75
          BZOUT=BZOUT*0.75
      ENDIF

      ELSE IF(XIN.LT.ULIMI(3)) THEN

      ENDIF

      b0halbasy=b0halbasyo
      xkhalbasy=xkhalbasyo
      zlhalbasy=zlhalbasyo
      zkhalbasy=zkhalbasyo
      fasym=fasymo
      ahwpol=ahwpolo
      xcenhal=xcenhalo

      RETURN
      END
