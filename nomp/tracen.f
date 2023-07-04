*CMZ :  2.61/02 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.57/00 25/10/2005  08.39.19  by  Michael Scheer
*CMZ :  2.54/05 01/06/2005  12.29.03  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE TRACEN(XCEN,YCEN,ZCEN,YPCEN,ZPCEN)
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

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,track0.
      include 'track0.cmn'
*KEND.

      DOUBLE PRECISION XCEN,YCEN,ZCEN,YPCEN,ZPCEN,
     &  XX,
     &  VZVX0,VZVX1,
     &  VYVX0,VYVX1,
     &  DX

      INTEGER I,I1,I2

      YCEN=-9999.
      YPCEN=-9999.
      ZCEN=-9999.
      ZPCEN=-9999.

      DX=WSXYZ(1,2)-WSXYZ(1,1)

      IF (XCEN.LT.WSXYZ(1,1)-DX.OR.XCEN.GT.WSXYZ(1,NCO)+DX) THEN
        RETURN
      ELSE IF (XCEN.LE.WSXYZ(1,1)) THEN
        I=1
      ELSE IF (XCEN.GE.WSXYZ(1,NCO)) THEN
        I=NCO-1
      ELSE
        I1=1
        I2=NCO
        DO WHILE (I2-I1.GT.1)
          I=(I2+I1)/2
          IF (XCEN.GE.WSXYZ(1,I)) THEN
            I1=I
          ELSE
            I2=I
          ENDIF
        ENDDO
        I=I1
      ENDIF

100   XX=(XCEN-WSXYZ(1,I))/(WSXYZ(1,I+1)-WSXYZ(1,I))

      YCEN=WSXYZ(2,I)+(WSXYZ(2,I+1)-WSXYZ(2,I))*XX

      VYVX0=WVXYZ(2,I)/WVXYZ(1,I)
      VYVX1=WVXYZ(2,I+1)/WVXYZ(1,I+1)

      YPCEN=VYVX0+(VYVX1-VYVX0)*XX

      ZCEN=WSXYZ(3,I)+(WSXYZ(3,I+1)-WSXYZ(3,I))*XX

      VZVX0=WVXYZ(3,I)/WVXYZ(1,I)
      VZVX1=WVXYZ(3,I+1)/WVXYZ(1,I+1)

      ZPCEN=VZVX0+(VZVX1-VZVX0)*XX

      RETURN
      END
