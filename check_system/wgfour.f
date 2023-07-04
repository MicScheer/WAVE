*CMZ :  3.02/08 25/06/2015  14.20.56  by  Michael Scheer
*CMZ :  3.02/05 23/03/2015  10.33.54  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.51/02 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.47/07 14/04/2003  12.01.12  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  16.30.02  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  2.14/02 27/04/2000  17.46.14  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.26.25  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  17.02.52  by  Michael Scheer
*CMZ : 00.02/00 11/12/96  14.30.51  by  Michael Scheer
*CMZ : 00.01/08 22/06/95  10.59.24  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.23.44  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.56  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.18  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WGFOUR
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
*KEND.

C--- CALCULATE FOURIER COEFFIENCE OF THE HORIZONTAL AND VERTICAL FOLDING
C    FUNCTIONS

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEND.

      INTEGER I,K,I1,IP,IM,NGCOEFP2,MFOUR,ISOUR,ICO
      INTEGER IGSIGZ,IGSIGY

      COMPLEX CKOEF(NGCOEFP/2+1+2)
      REAL*4  YFOUR(NGCOEFP+2+2),AKOEF(NGCOEFP/2+1+2)
      EQUIVALENCE (CKOEF,YFOUR)

      DOUBLE PRECISION XFOUR(NGCOEFP+2)
      DOUBLE PRECISION XLFOUR,DXFOUR,BY

      DOUBLE PRECISION FOUFUNX
      DOUBLE PRECISION DGSIGZO,DGSIGYO

C--- LOOP UBER NGCOEFP-PUNKTE FUER FAST-FOURIER-TRANSFORMATION

C--- DETERMINE INTERVALL FOR FFT I.E. COMPARE SPACING OF GRID
C    AND SIGMAS OF FOLDING FUNCTION AND ADJUST DGSIG

      IF (NGFOURZ.GT.NGCOEFP/2) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN WGFOUR:'
          WRITE(LUNGFO,*)'NGFOURZ.GT.NGCOEFP/2'
          WRITE(LUNGFO,*)
     &'CHECK NGFOURZ IN NAMELIST WFOLDN OR INCREASE NGCOEFP IN CMPARA.CMN'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN WGFOUR:'
          WRITE(6,*)'NGFOURZ.GT.NGCOEFP/2'
          WRITE(6,*)
     &'CHECK NGFOURZ IN NAMELIST WFOLDN OR INCREASE NGCOEFP IN CMPARA.CMN'
          WRITE(6,*)
          STOP '*** PROGRAMM ABORTED  ***'
      ENDIF

      IF (NGFOURY.GT.NGCOEFP/2) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN WGFOUR:'
          WRITE(LUNGFO,*)'NGFOURY.GT.NGCOEFP/2'
          WRITE(LUNGFO,*)
     &'CHECK NGFOURY IN NAMELIST WFOLDN OR INCREASE NGCOEFP IN CMPARA.CMN'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN WGFOUR:'
          WRITE(6,*)'NGFOURY.GT.NGCOEFP/2'
          WRITE(6,*)
     &'CHECK NGFOURY IN NAMELIST WFOLDN OR INCREASE NGCOEFP IN CMPARA.CMN'
          WRITE(6,*)
          STOP '*** PROGRAMM ABORTED  ***'
      ENDIF

C--- LOOP OVER ALL SOURCES

      DO ISOUR=1,NSOURCE

      IF (IF1DIM.EQ.0.AND.WSIGZ(ISOUR).EQ.0.D0) THEN
          WRITE(6,*)'*** ERROR IN WGFOUR ***'
          WRITE(6,*)'*** SIGMA FOR HORIZONTALFOLDING FUNCTION IS ZERO ***'
          WRITE(6,*)'*** SOURCE:',ISOUR
          WRITE(6,*)'*** TOTAL NUMBER OF SOURCES:',NSOURCE
          WRITE(6,*)'*** CHECK NAMELIST WFOLDN IN INPUT FILE ***'
          WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
          WRITE(LUNGFO,*)'*** ERROR IN WGFOUR ***'
          WRITE(LUNGFO,*)'*** SIGMA FOR HORIZONTAL FOLDING FUNCTION IS ZERO ***'
          WRITE(LUNGFO,*)'*** SOURCE:',ISOUR
          WRITE(LUNGFO,*)'*** TOTAL NUMBER OF SOURCES:',NSOURCE
          WRITE(LUNGFO,*)'*** CHECK NAMELIST WFOLDN IN INPUT FILE ***'
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
          STOP
      ENDIF !(WSIGZ(ISOUR).EQ.0.D0)

      IF (WSIGY(ISOUR).EQ.0.D0) THEN
          WRITE(6,*)'*** ERROR IN WGFOUR ***'
          WRITE(6,*)'*** SIGMA FOR VERTICAL FOLDING FUNCTION IS ZERO ***'
          WRITE(6,*)'*** SOURCE:',ISOUR
          WRITE(6,*)'*** TOTAL NUMBER OF SOURCES:',NSOURCE
          WRITE(6,*)'*** CHECK NAMELIST WFOLDN IN INPUT FILE ***'
          WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
          WRITE(LUNGFO,*)'*** ERROR IN WGFOUR ***'
          WRITE(LUNGFO,*)'*** SIGMA FOR VERTICAL FOLDING FUNCTION IS ZERO ***'
          WRITE(LUNGFO,*)'*** SOURCE:',ISOUR
          WRITE(LUNGFO,*)'*** TOTAL NUMBER OF SOURCES:',NSOURCE
          WRITE(LUNGFO,*)'*** CHECK NAMELIST WFOLDN IN INPUT FILE ***'
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
          STOP
      ENDIF !(WSIGZ(ISOUR).EQ.0.D0)

      IGSIGZ=NINT(DSIGZ(ISOUR)/OBSVDZ)
      IGSIGY=NINT(DSIGY(ISOUR)/OBSVDY)

C        HALF WIDTH OF FUNCTION ON BOTTOM LINE
      DSIGZ(ISOUR)=DFLOAT(IGSIGZ+1)*OBSVDZ
      DSIGY(ISOUR)=DFLOAT(IGSIGY+1)*OBSVDY

      DSIG2Z(ISOUR)=DSIGZ(ISOUR)**2
      DSIG2Y(ISOUR)=DSIGY(ISOUR)**2

      DGSIGZO=DGSIGZ(ISOUR)
      DGSIGYO=DGSIGY(ISOUR)

      IF (IF1DIM.EQ.0) DGSIGZ(ISOUR)=DSIGZ(ISOUR)/WSIGZ(ISOUR)
      DGSIGY(ISOUR)=DSIGY(ISOUR)/WSIGY(ISOUR)

      IF (DGSIGZ(ISOUR)-DGSIGZO.GT.1.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING IN WGFOUR:'
          WRITE(LUNGFO,*)
     &'ADJUSTED VALUES OF DGSIGZ(ISOUR) DIFFER MORE THEN SIGMA FROM GIVEN VALUE'
          WRITE(LUNGFO,*)'CHECK FOLDING FUNCTION AND DGSIGZ'
c20150625          WRITE(LUNGFO,*)'DGSIGZ LIMITED TO 10*WSIGZ'
c          WRITE(LUNGFO,*)
c          WRITE(6,*)
c          WRITE(6,*)'*** WARNING IN WGFOUR:'
c          WRITE(6,*)
c     &'ADJUSTED VALUES OF DGSIGZ(ISOUR) DIFFER MORE THEN SIGMA FROM GIVEN VALUE'
c          WRITE(6,*)'CHECK FOLDING FUNCTION AND DGSIGZ'
c          WRITE(6,*)
c20150625          dgsigz(isour)=10.0d0*wsigz(isour)
      ENDIF

      IF (DGSIGY(ISOUR)-DGSIGYO.GT.1.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR WGFOUR:'
          WRITE(LUNGFO,*)
     &'ADJUSTED VALUES OF DGSIGY(ISOUR) DIFFER MORE THEN SIGMA FROM GIVEN VALUE'
          WRITE(LUNGFO,*)'CHECK FOLDING FUNCTION AND DGSIGY'
          WRITE(LUNGFO,*)
c20150625          WRITE(LUNGFO,*)'DGSIGY LIMITED TO 10*WSIGY'
c          WRITE(6,*)
c          WRITE(6,*)'*** WARNING SR WGFOUR:'
c          WRITE(6,*)
c     &'ADJUSTED VALUES OF DGSIGY(ISOUR) DIFFER MORE THEN SIGMA FROM GIVEN VALUE'
c          WRITE(6,*)'CHECK FOLDING FUNCTION AND DGSIGY'
c          WRITE(6,*)
c20150625          dgsigy(isour)=10.0d0*wsigy(isour)
      ENDIF

C--- VERTICAL DIRECTION

      XLFOUR=2.*DSIGY(ISOUR)

      DXFOUR=XLFOUR/NGCOEFP
      NGCOEFP2=NGCOEFP/2
      MFOUR=NINT(ALOG(FLOAT(NGCOEFP))/ALOG(2.E0))

      DO I=1,NGCOEFP2+1          !SYMMETRISCHE X-WERTE
          XFOUR(I)          =-DXFOUR*(NGCOEFP2+1-I)
          XFOUR(NGCOEFP+1-I+1)=-XFOUR(I)
      END DO

      DO I=1,NGCOEFP2+1          !SYMMETRISCHE Y-WERTE

          I1=I-1
          IP=NGCOEFP2+1+I1
          IM=NGCOEFP2+1-I1

          BY=FOUFUNX(XFOUR(IP),WSIGY(ISOUR))

          YFOUR(IP)=BY
          YFOUR(IM)=BY

      END DO

      CALL RFFT(CKOEF,-MFOUR) !FFT MIT CERN-ROUTINE D703

      DO K=1,NGCOEFP2+1 !REELLE KOEFFIZIENTEN
          AKOEF(K)=(-1.)**(K-1)*2.*REAL(CKOEF(K))
      ENDDO

C--- CALCULATE COEFFICIENTS FOR FOLDING FUNCTION

      DO ICO=1,NGFOURY

          YKGAUSS(ICO,ISOUR)=2.D0*PI1/XLFOUR*DFLOAT(ICO)
          GCOEFV(ICO,ISOUR)=AKOEF(ICO)
          IF (ICO.EQ.1) GCOEFV(ICO,ISOUR)=GCOEFV(ICO,ISOUR)*0.5D0

      ENDDO !NGFOURY

C--- HORIZONTAL DIRECTION

      IF (IF1DIM.EQ.0) THEN

      XLFOUR=2.*DSIGZ(ISOUR)

      DXFOUR=XLFOUR/NGCOEFP
      NGCOEFP2=NGCOEFP/2
      MFOUR=NINT(ALOG(FLOAT(NGCOEFP))/ALOG(2.E0))

      DO I=1,NGCOEFP2+1          !SYMMETRISCHE X-WERTE
          XFOUR(I)          =-DXFOUR*(NGCOEFP2+1-I)
          XFOUR(NGCOEFP+1-I+1)=-XFOUR(I)
      END DO

      DO I=1,NGCOEFP2+1          !SYMMETRISCHE Y-WERTE

          I1=I-1
          IP=NGCOEFP2+1+I1
          IM=NGCOEFP2+1-I1

          BY=FOUFUNX(XFOUR(IP),WSIGZ(ISOUR))

          YFOUR(IP)=BY
          YFOUR(IM)=BY

      END DO

      CALL RFFT(CKOEF,-MFOUR) !FFT MIT CERN-ROUTINE D703

      DO K=1,NGCOEFP2+1 !REELLE KOEFFIZIENTEN
          AKOEF(K)=(-1.)**(K-1)*2.*REAL(CKOEF(K))
      ENDDO

C--- CALCULATE COEFFICIENTS FOR FOLDING FUNCTION

      DO ICO=1,NGFOURZ

          XKGAUSS(ICO,ISOUR)=2.D0*PI1/XLFOUR*DFLOAT(ICO)
          GCOEFH(ICO,ISOUR)=AKOEF(ICO)
          IF (ICO.EQ.1) GCOEFH(ICO,ISOUR)=GCOEFH(ICO,ISOUR)*0.5D0
      ENDDO !NGFOURZ
      ENDIF !IF1DIM

      ENDDO !ISOUR

      RETURN
      END
