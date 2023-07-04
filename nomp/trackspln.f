*CMZ :  2.63/03 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.16/08 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  16.15.52  by  Michael Scheer
*CMZ :  2.10/01 17/02/99  15.53.29  by  Michael Scheer
*-- Author :    Michael Scheer   16/02/99

      SUBROUTINE TRACKSPLN(X1,DTIM,
     &                  X2,Y2,Z2,VX2,VY2,VZ2)
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

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,track.
      include 'track.cmn'
*KEND.

      INTEGER ICAL,I,ICOUNT

      DOUBLE PRECISION X1,DTIM,X2,Y2,Z2,VX2,VY2,VZ2,T1,T2,XDUM

      DOUBLE PRECISION S2T(NWMAXP),S2X(NWMAXP),S2Y(NWMAXP),S2Z(NWMAXP)
      DOUBLE PRECISION S2VX(NWMAXP),S2VY(NWMAXP),S2VZ(NWMAXP)
      DOUBLE PRECISION WS1(NWMAXP),WS2(NWMAXP),WS3(NWMAXP),WS4(NWMAXP)

      DOUBLE PRECISION XA(NWMAXP),YA(NWMAXP),ZA(NWMAXP)
      DOUBLE PRECISION VXA(NWMAXP),VYA(NWMAXP),VZA(NWMAXP)

      DATA ICAL/0/

      WRITE(6,*)'*** ERROR IN TRACKSPLN ***'
      WRITE(6,*)'*** ROUTINE NOT YET ADAPTED TO F90 ***'
      STOP

      IF (ICAL.EQ.0) THEN

            DO I=1,NCO

         XA(I)=WSXYZ(1,I)
         YA(I)=WSXYZ(2,I)
         ZA(I)=WSXYZ(3,I)

         VXA(I)=WVXYZ(1,I)
         VYA(I)=WVXYZ(2,I)
         VZA(I)=WVXYZ(3,I)

            ENDDO

            CALL util_spline_coef(XA,WTIM0,NCO,-9999.0d0,-9999.0d0,S2T,WS1,WS2,WS3,WS4)
          CALL TIME_SPLINE_INTER(XA,WTIM0,S2T,NCO,X1,T1,-1,ICOUNT)

            CALL util_spline_coef(WTIM0,XA,NCO,-9999.0d0,-9999.0d0,S2X,WS1,WS2,WS3,WS4)
            CALL util_spline_coef(WTIM0,YA,NCO,-9999.0d0,-9999.0d0,S2Y,WS1,WS2,WS3,WS4)
            CALL util_spline_coef(WTIM0,ZA,NCO,-9999.0d0,-9999.0d0,S2Z,WS1,WS2,WS3,WS4)

            CALL util_spline_coef(WTIM0,VXA,NCO,-9999.0d0,-9999.0d0,S2VX,WS1,WS2,WS3,WS4)
            CALL util_spline_coef(WTIM0,VYA,NCO,-9999.0d0,-9999.0d0,S2VY,WS1,WS2,WS3,WS4)
            CALL util_spline_coef(WTIM0,VZA,NCO,-9999.0d0,-9999.0d0,S2VZ,WS1,WS2,WS3,WS4)

          CALL XYZ_SPLINE_INTER(WTIM0,XA,S2X,NCO,T1,XDUM,-1,ICOUNT)

          ICAL=1

      ENDIF !ICAL

C--- GET T1 OF TIME X1

      CALL TIME_SPLINE_INTER(XA,WTIM0,S2T,NCO,X1,T1,0,ICOUNT)

C--- GET X2,Y2,Z2 OF TIME T2

      T2=T1+DTIM

      CALL XYZ_SPLINE_INTER(WTIM0,XA,S2X,NCO,T2,X2,0,ICOUNT)
      CALL XYZ_SPLINE_INTER(WTIM0,YA,S2Y,NCO,T2,Y2,0,ICOUNT)
      CALL XYZ_SPLINE_INTER(WTIM0,ZA,S2Z,NCO,T2,Z2,0,ICOUNT)

      CALL XYZ_SPLINE_INTER(WTIM0,VXA,S2VX,NCO,T2,VX2,0,ICOUNT)
      CALL XYZ_SPLINE_INTER(WTIM0,VYA,S2VY,NCO,T2,VY2,0,ICOUNT)
      CALL XYZ_SPLINE_INTER(WTIM0,VZA,S2VZ,NCO,T2,VZ2,0,ICOUNT)

      ICOUNT=0

      RETURN
      END



