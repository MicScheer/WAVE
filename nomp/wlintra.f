*CMZ :  3.00/00 02/04/2013  13.39.05  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/09 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.66/07 17/12/2009  13.18.30  by  Michael Scheer
*CMZ :  2.16/08 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.57.05  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.12  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WLINTRA
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
*KEEP,wbetaf90u.
      include 'wbetaf90u.cmn'
*KEND.

C--- CALCULATES LINEARE TRANSFER MATRICES WITHOUT COUPLING FROM BETA-FUNCTION
C    STORED IN ARRAY WBETA AND TUNE SHIFTS

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,wbetaf90.
      include 'wbetaf90.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER IP,ical

      DOUBLE PRECISION ALPHA0H,ALPHAH,BETA0H,BETAH,SINPHIH,COSPHIH
      DOUBLE PRECISION ALPHA0V,ALPHAV,BETA0V,BETAV,SINPHIV,COSPHIV
      DOUBLE PRECISION T11,T12,BY,BETA

      data ical/0/

      if (ical.ne.0) return
      ical=1

C--- TUNE SHIFT

      TUNEH=0.0
      TUNEV=0.0

      WTUNE(1,1)=TUNEH
      WTUNE(2,1)=TUNEV

      DO IP=2,NCO

          TUNEH=TUNEH+2.D0/(WBETA(2,IP)+WBETA(2,IP-1))*DS0
          TUNEV=TUNEV+2.D0/(WBETA(4,IP)+WBETA(4,IP-1))*DS0

          WTUNE(1,IP)=TUNEH
          WTUNE(2,IP)=TUNEV

      ENDDO !IP


C--- LINEAR TRANSFERMATRICES

      BETA0H=WBETA(2,1)
      BETA0V=WBETA(4,1)
      ALPHA0H=-WBETA(3,1)/2.D0
      ALPHA0V=-WBETA(5,1)/2.D0

      DO IP=1,NCO

         TUNEH=WTUNE(1,IP)
         TUNEV=WTUNE(2,IP)

         COSPHIH=DCOS(TUNEH)
         COSPHIV=DCOS(TUNEV)
         SINPHIH=DSIN(TUNEH)
         SINPHIV=DSIN(TUNEV)
         BETAH  =WBETA(2,IP)
         BETAV  =WBETA(4 ,IP)
         ALPHAH =-WBETA(3,IP)/2.D0
         ALPHAV =-WBETA(5,IP)/2.D0

         WLTM(1,1,IP)=DSQRT(BETAH/BETA0H)*(COSPHIH+ALPHA0H*SINPHIH)
         WLTM(1,2,IP)=DSQRT(BETA0H*BETAH)*SINPHIH
         WLTM(1,3,IP)=1.0D0/(DSQRT(BETA0H*BETAH))
     &     *(
     &     (ALPHA0H-ALPHAH) *COSPHIH
     &     -
     &     (1.D0+ALPHA0H*ALPHAH) *SINPHIH
     &     )
         WLTM(1,4,IP)=DSQRT(BETA0H/BETAH)*(COSPHIH-ALPHAH*SINPHIH)

         WLTM(2,1,IP)=DSQRT(BETAV/BETA0V)*(COSPHIV+ALPHA0V*SINPHIV)
         WLTM(2,2,IP)=DSQRT(BETA0V*BETAV)*SINPHIV
         WLTM(2,3,IP)=1.0D0/(DSQRT(BETA0V*BETAV))
     &     *(
     &     (ALPHA0V-ALPHAV) *COSPHIV
     &     -
     &     (1.D0+ALPHA0V*ALPHAV) *SINPHIV
     &     )
         WLTM(2,4,IP)=DSQRT(BETA0V/BETAV)*(COSPHIV-ALPHAV*SINPHIV)

      ENDDO !IP

      TMH(1,1)=WLTM(1,1,NCO)
      TMH(1,2)=WLTM(1,2,NCO)
      TMH(2,1)=WLTM(1,3,NCO)
      TMH(2,2)=WLTM(1,4,NCO)

      TMV(1,1)=WLTM(2,1,NCO)
      TMV(1,2)=WLTM(2,2,NCO)
      TMV(2,1)=WLTM(2,3,NCO)
      TMV(2,2)=WLTM(2,4,NCO)

C--- VERTICAL TUNE SHIFT

      TUNSHI=0.0
      DO IP=1,NCO
        BY=WTRA(2,3,IP)
        BETA=WBETA(4,IP)
        TUNSHI=TUNSHI+BETA*BY**2
      ENDDO !IP
      TUNSHI=TUNSHI*DS0/(4.*PI1)/DMYENERGY**2*0.3**2

C--- TUNES (PHASE ADVANCE) OF CORRESPONDING DRIFT

      T11=1.D0
      T12=WTRA(1,1,NCO)-WTRA(1,1,1)
      TUNEH0=DATAN(T12/(T11*WBETA(2,1)+WBETA(3,1)/2.*T12))
      TUNEV0=DATAN(T12/(T11*WBETA(4,1)+WBETA(5,1)/2.*T12))
      IF(TUNEH0.LT.0.0)TUNEH0=TUNEH0+PI1
      IF(TUNEV0.LT.0.0)TUNEV0=TUNEV0+PI1

      RETURN
      END
