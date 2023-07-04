*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.56/00 25/06/2010  12.15.46  by  Michael Scheer
*CMZ :  2.55/00 10/08/2005  15.45.45  by  Michael Scheer
*CMZ :  2.51/00 14/05/2004  14.43.12  by  Michael Scheer
*CMZ :  2.47/20 01/12/2003  15.25.47  by  Michael Scheer
*CMZ :  2.37/04 05/12/2001  17.36.28  by  Michael Scheer
*CMZ :  2.34/09 18/09/2001  22.01.04  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  16.15.48  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.16/00 09/06/2000  12.29.25  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  2.13/10 14/04/2000  14.31.40  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.43.05  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.17  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WSIGFOL
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

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEEP,wbetaf90u.
      include 'wbetaf90u.cmn'
*KEND.

C--- COMPUTES RMS OF 2D GAUSSIAN FOR FOLDING PROCEDURE
C    THE SIGMAS DEPEND ON DISTANCE OF OBSERVER AND SOURCE POINT

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,wbetaf90.
      include 'wbetaf90.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEND.

      INTEGER ISOUR,ICEN

      DOUBLE PRECISION BETH,BETPH,BETV,BETPV,DIST2,DIST,ALFH,ALFV,GAMH,GAMV

      IF(NSOURCE.GT.LIDIMP) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)' *** ERROR IN WSIGFOL ***'
        WRITE(LUNGFO,*)'NSOURCE.GT.LIDIMP'
        WRITE(LUNGFO,*)'CHECK INPUT FILE WAVE.IN OR INCREASE DIMENSION'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)' *** ERROR IN WSIGFOL ***'
        WRITE(6,*)'NSOURCE.GT.LIDIMP'
        WRITE(6,*)'CHECK INPUT FILE WAVE.IN OR INCREASE DIMENSION'
        WRITE(6,*)
        STOP
      ENDIF

C--- CENTER OF SOURCE POINTS, BETA-FUNCTION

      DO ISOUR=1,NSOURCE

        IF (BSIGZ(ISOUR).NE.0.D0.OR.BSIGZP(ISOUR).NE.0.D0) THEN
          WSIGZ(ISOUR)=SQRT(BSIGZ(ISOUR)**2+(PINCEN(1)*BSIGZP(ISOUR))**2)
        ENDIF

        IF (BSIGY(ISOUR).NE.0.D0.OR.BSIGYP(ISOUR).NE.0.D0) THEN
          WSIGY(ISOUR)=SQRT(BSIGY(ISOUR)**2+(PINCEN(1)*BSIGYP(ISOUR))**2)
        ENDIF

        IF (WSIGZ(ISOUR).EQ.0.D0) THEN
          WRITE(LUNGFO,*)
     &      '*** WARNING IN WSIGFOL: WSIGZ(',ISOUR,') taken from WSIGZ(1)'
          WRITE(LUNGFO,*)
     &     '*** WILL BE WRONG IF DISTANCE FROM SOURCE TO OBSERVER DIFFERS'
          WSIGZ(ISOUR)=WSIGZ(1)
        ENDIF

        IF (WSIGY(ISOUR).EQ.0.D0) THEN
          WRITE(LUNGFO,*)
     &      '*** WARNING IN WSIGFOL: WSIGY(',ISOUR,') taken from WSIGY(1)'
          WRITE(LUNGFO,*)
     &     '*** WILL BE WRONG IF DISTANCE FROM SOURCE TO OBSERVER DIFFERS'
          WSIGY(ISOUR)=WSIGY(1)
        ENDIF

C--- USER DEFINED SIGMAS

        IF (ISIGUSR.NE.0) THEN

          WSIG2Z(ISOUR)=WSIGZ(ISOUR)**2
          WSIG2Y(ISOUR)=WSIGY(ISOUR)**2

        ELSE

C- DISTANCE FROM PINHOLE CENTER TO CENTER OF SOURCE

          DIST2=(SOURCEN(1,1,ISOUR)-PINCEN(1))**2
     &      +(SOURCEN(2,1,ISOUR)-PINCEN(2))**2
     &      +(SOURCEN(3,1,ISOUR)-PINCEN(3))**2
          DIST=DSQRT(DIST2)

C- BETA-FUNCTION AND DERIVATIVE

          ICEN=ISOURCEN(ISOUR)
          BETH=WBETA(2,ICEN)
          BETPH=WBETA(3,ICEN)
          BETV=WBETA(4,ICEN)
          BETPV=WBETA(5,ICEN)
          ALFH=-BETPH/2.
          GAMH=(1.+ALFH**2)/BETH
          ALFV=-BETPV/2.
          GAMV=(1.+ALFV**2)/BETV

          WSIG2Z(ISOUR)=EPS0H*(BETH+DIST2*GAMH-2.*DIST*ALFH)
          WSIG2Y(ISOUR)=EPS0V*(BETV+DIST2*GAMV-2.*DIST*ALFV)


          WSIGZ(ISOUR)=DSQRT(WSIG2Z(ISOUR))
          WSIGY(ISOUR)=DSQRT(WSIG2Y(ISOUR))

        ENDIF !ISGUSR

        IF (IF1DIM.EQ.1) THEN
          WSIGZ(ISOUR)=0.0
          DGSIGZ(ISOUR)=0.0
        ENDIF  !IF1DIM

        IF (DGSIGZ(ISOUR).EQ.0.D0) DGSIGZ(ISOUR)=DGSIGZ(1)
        IF (DGSIGY(ISOUR).EQ.0.D0) DGSIGY(ISOUR)=DGSIGY(1)

        DSIGZ(ISOUR)=WSIGZ(ISOUR)*DGSIGZ(ISOUR)
        DSIGY(ISOUR)=WSIGY(ISOUR)*DGSIGY(ISOUR)

        IF (  IF1DIM.EQ.0 .AND.
     &      (WSIGZ(ISOUR).EQ.0.0
     &      .OR.  DGSIGZ(ISOUR).EQ.0.0)
     &      ) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)' *** ERROR IN WSIGFOL ***'
          WRITE(LUNGFO,*)'ZERO WIDTH OF HORIZONTAL FOLDING FUNCTION'
          WRITE(LUNGFO,*)'CHECK WSIGZ AND DGSIGZ IN NAMELIST WFOLDN'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)' *** ERROR IN WSIGFOL ***'
          WRITE(6,*)'ZERO WIDTH OF HORIZONTAL FOLDING FUNCTION'
          WRITE(6,*)'CHECK WSIGZ AND DGSIGZ IN NAMELIST WFOLDN'
          WRITE(6,*)
          STOP
        ENDIF

        IF (
     &      WSIGY(ISOUR).EQ.0.0
     &      .OR.  DGSIGY(ISOUR).EQ.0.0
     &      ) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)' *** ERROR IN WSIGFOL ***'
          WRITE(LUNGFO,*)'ZERO WIDTH OF VERTICAL FOLDING FUNCTION'
          WRITE(LUNGFO,*)'CHECK WSIGY AND DGSIGY IN NAMELIST WFOLDN'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)' *** ERROR IN WSIGFOL ***'
          WRITE(6,*)'ZERO WIDTH OF VERTICAL FOLDING FUNCTION'
          WRITE(6,*)'CHECK WSIGY AND DGSIGY IN NAMELIST WFOLDN'
          WRITE(6,*)
          STOP
        ENDIF

      ENDDO !ISOUR

      IF(ISTOKES.NE.0 .AND. ISIGSTO.GT.NSOURCE) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)' *** ERROR IN WSIGFOL ***'
        WRITE(LUNGFO,*)'ISIGSTO GREATER THAN NUMBER OF SOURCE POINTS'
        WRITE(LUNGFO,*)'CHECK INPUT FILE WAVE.IN'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)' *** ERROR IN WSIGFOL ***'
        WRITE(6,*)'ISIGSTO GREATER THAN NUMBER OF SOURCE POINTS'
        WRITE(6,*)'CHECK INPUT FILE WAVE.IN'
        WRITE(6,*)
        STOP
      ENDIF

      RETURN
      END
