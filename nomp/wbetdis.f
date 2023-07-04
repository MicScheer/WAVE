*CMZ :  3.04/00 11/01/2018  10.11.29  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.09.17  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.66/13 07/07/2010  11.05.47  by  Michael Scheer
*CMZ :  2.66/09 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.66/07 20/01/2010  16.23.37  by  Michael Scheer
*CMZ :  2.47/12 17/12/2009  14.02.55  by  Michael Scheer
*CMZ :  2.41/13 22/08/2002  17.20.15  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  16.19.29  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.43.04  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  10.41.10  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.04  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.13  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WBETDIS
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

C--- CALCULATE BETA-FUNCTIONS, DISPERSION AND THEIR DERIVATIVES
C    WRITES RESULTS TO FILES

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,wbetaf90.
      include 'wbetaf90.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER IP,IW,iwarn
      DOUBLE PRECISION DIS,DISP,HBET,HBETP,VBET,VBETP

      data iwarn/0/

      ALLOCATE(WBETA(16,NCO))
      ALLOCATE(WBZZPYYP(4,NCO))
      ALLOCATE(WBETAK(3,NCO))
      ALLOCATE(WLTM(2,4,NCO))
      ALLOCATE(WTUNE(2,NCO))
      wbeta=0.0d0
      wbzzpyyp=0.0d0
      wbetak=0.0d0
      wltm=0.0d0
      wtune=0.0d0

C--- CHECK PLANARITY OF REFERENCE ORBIT

      DO IP=1,NCO

        IF (abs(WSXYZ(2,IP)).gt.1.0d-10.and.iwarn.eq.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING IN WBETDIS ***'
          WRITE(LUNGFO,*)'REFERENCE ORBIT NOT PLANAR'
          WRITE(LUNGFO,*)
          WRITE(6,*) '*** WARNING IN WBETDIS ***'
          WRITE(6,*)'REFERENCE ORBIT NOT PLANAR'
          WRITE(6,*)
          iwarn=1
c          STOP
        ENDIF

      ENDDO

C--- BETA FUNCTION

      CALL WBETFN

C--- LINEARE TRANSFER MATRICES

      CALL WLINTRA

C--- DISPERSION

      CALL WDISPER

C--- RADIATION INTEGRALS

      CALL WI2I4I5

C--- TRAJECTORY IN PHASE SPACE

      CALL WZZPYYP

C--- WRITE RESULTS TO FILE


C    WBETA(2,*) HORIZONTAL BETA-FUNCTION
C    WBETA(3,*) DERIVATIVE OF HORIZONTAL BETA-FUNCTION
C    WBETA(4,*) VERTICAL BETA-FUNCTION
C    WBETA(5,*) DERIVATIVE OF VERTICAL BETA-FUNCTION

C    WBETA(6,*) DISPERTION (SR WDISPER)
C    WBETA(7,*) DERIVATIVE OF DISPERTION (SR WDISPER)

C    OPTICAL FUNCTIONS ALPA, BETA AND GAMMA
C         ALPHA=-WBETA(1+IT*2,IP)/2.D0
C         !TEMPORARY VALUE OF WBETAK(1/2,IP) IS DERIVATIVE OF WBETA(2/4,IP)
C         ALPHAP=-WBETAK(IT,IP)/2.D0
C         BETA=WBETA(2*IT,IP)
C         GAMA=(1+ALPHA**2)/BETA
C         WBETAK(IT,IP)=(ALPHAP+GAMA)/BETA
C--- RECALCULATE MAGNETIC FIELD FROM Ky+Kz=1/rho**2
C         WBETAK(3,IP)=DSQRT(WBETAK(2,IP)+WBETAK(1,IP))*EMOM/CLIGHT1

      OPEN(UNIT=LUNWB,FILE=FILEWB,STATUS='unknown',recl=256)

      DO IP=1,NCO
        WRITE(LUNWB,*)ICODE,(SNGL(WBETA(IW,IP)),IW=1,9),
     &    (SNGL(WBETAK(IW,IP)),IW=1,3)
      ENDDO   !NCO

      CLOSE(LUNWB)

C--- HBOOK

      IF (IHBETA.NE.0) CALL HBETA

C--- FUNCTIONS AT CENTER OF WLS

         HBET=WBETA(2,NCO/2)
         HBETP=WBETA(3,NCO/2)
         VBET=WBETA(4,NCO/2)
         VBETP=WBETA(5,NCO/2)
         DIS=WBETA(6,NCO/2)
         DISP=WBETA(7,NCO/2)

C090792  DO IP=1,NCO-1
C
C         IF (WBETA(1,IP)*WBETA(1,IP+1).LE.0.0) THEN
C        HBET=WBETA(2,IP)
C     &         +(WBETA(2,IP+1)-WBETA(2,IP))
C     &         /(WBETA(1,IP+1)-WBETA(1,IP))
C     &         *(0.0-WBETA(1,IP))
C        HBETP=WBETA(3,IP)
C     &         +(WBETA(3,IP+1)-WBETA(3,IP))
C     &         /(WBETA(1,IP+1)-WBETA(1,IP))
C     &         *(0.0-WBETA(1,IP))
C        VBET=WBETA(4,IP)
C     &         +(WBETA(4,IP+1)-WBETA(4,IP))
C     &         /(WBETA(1,IP+1)-WBETA(1,IP))
C     &         *(0.0-WBETA(1,IP))
C        VBETP=WBETA(5,IP)
C     &         +(WBETA(5,IP+1)-WBETA(5,IP))
C     &         /(WBETA(1,IP+1)-WBETA(1,IP))
C     &         *(0.0-WBETA(1,IP))
C        DIS=WBETA(6,IP)
C     &         +(WBETA(6,IP+1)-WBETA(6,IP))
C     &         /(WBETA(1,IP+1)-WBETA(1,IP))
C     &         *(0.0-WBETA(1,IP))
C        DISP=WBETA(7,IP)
C     &         +(WBETA(7,IP+1)-WBETA(7,IP))
C     &         /(WBETA(1,IP+1)-WBETA(1,IP))
C     &         *(0.0-WBETA(1,IP))
C        GOTO 90
C         ENDIF
C
C     ENDDO !IP


c90    WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  'HOR. BETA-FUNCTION AND DERIVATIVE AT ENTRANCE OF WLS'
      WRITE(LUNGFO,*)SNGL(WBETA(2,1)),SNGL(WBETA(3,1))
      WRITE(LUNGFO,*)
     &  'VER. BETA-FUNCTION AND DERIVATIVE AT ENTRANCE OF WLS'
      WRITE(LUNGFO,*)SNGL(WBETA(4,1)),SNGL(WBETA(5,1))
      WRITE(LUNGFO,*)
     &  'HOR. BETA-FUNCTION AND DERIVATIVE AT CENTER OF WLS'
      WRITE(LUNGFO,*)SNGL(HBET),SNGL(HBETP)
      WRITE(LUNGFO,*)
     &  'VER. BETA-FUNCTION AND DERIVATIVE AT CENTER OF WLS'
      WRITE(LUNGFO,*)SNGL(VBET),SNGL(VBETP)
      WRITE(LUNGFO,*)
     &  'DISPERSION AND DERIVATIVE AT CENTER OF WLS'
      WRITE(LUNGFO,*)SNGL(DIS),SNGL(DISP)
      WRITE(LUNGFO,*)
     &  'HOR. BETA-FUNCTION AND DERIVATIVE AT EXIT OF WLS'
      WRITE(LUNGFO,*)SNGL(WBETA(2,NCO)),SNGL(WBETA(3,NCO))
      WRITE(LUNGFO,*)
     &  'VER. BETA-FUNCTION AND DERIVATIVE AT EXIT OF WLS'
      WRITE(LUNGFO,*)SNGL(WBETA(4,NCO)),SNGL(WBETA(5,NCO))
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  'TOTAL HORIZONTAL AND VERTICAL PHASE ADVANCES:'
      WRITE(LUNGFO,*)SNGL(TUNEH),SNGL(TUNEV)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  'TOTAL HORIZONTAL AND VERTICAL PHASE ADVANCES FOR CORRESPONDING DRIFT:'
      WRITE(LUNGFO,*)SNGL(TUNEH0),SNGL(TUNEV0),
     &                 '(SAME BETA AT STARTING POINT)'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  'VERTICAL TUNE SHIFT (FROM 1/4/PI*INTEGRAL(BETA/RHO**2)):',
     &   SNGL(TUNSHI)
      WRITE(LUNGFO,*)
     &   '(ONLY CORRECT FOR INSERTION DEVICES WITHOUT TRANSVERSAL GRADIENT,'
      WRITE(LUNGFO,*)
     &   'NOT CORRECT E.G. FOR QUADRUPOLS)'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  'HORIZONTAL LINEAR TRANSFERMATRIX CALCULATED FROM BETA-FUNCTION:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)SNGL(TMH(1,1)),SNGL(TMH(1,2))
      WRITE(LUNGFO,*)SNGL(TMH(2,1)),SNGL(TMH(2,2))
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  'VERTICAL LINEAR TRANSFERMATRIX CALCULATED FROM BETA-FUNCTION:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)SNGL(TMV(1,1)),SNGL(TMV(1,2))
      WRITE(LUNGFO,*)SNGL(TMV(2,1)),SNGL(TMV(2,2))
      WRITE(LUNGFO,*)

      RETURN
      END
