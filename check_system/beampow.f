*CMZ :  4.00/14 30/12/2021  15.41.22  by  Michael Scheer
*CMZ :  4.00/13 07/12/2021  18.47.10  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.15.16  by  Michael Scheer
*CMZ :  3.02/04 23/01/2015  16.14.38  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  10.40.59  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  11.18.26  by  Michael Scheer
*CMZ :  2.67/01 13/03/2012  12.31.13  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  15.15.38  by  Michael Scheer
*CMZ :  2.56/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.51/00 24/05/2004  15.02.20  by  Michael Scheer
*CMZ :  2.35/02 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.16/08 25/10/2000  16.44.25  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  1.03/01 10/02/98  17.32.59  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.27  by  Michael Scheer
*CMZ : 00.01/06 13/02/95  10.20.34  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  14.13.58  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.46.53  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.26  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BEAMPOW(NPOL)

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

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

C--- CALCULATES POWER DISTRIBUTION ALONG THE WALLS OF THE BEAMLINE
C    STORES RESULTS IN RADPOW(19,NDPOLP,NWMAXP)
C
C     RADPOW(1:2,..)   X-COORDINATE 1. AND 2. WALL
C     RADPOW(3:4,..)   POWER DENSITY
C     RADPOW(5:6,..)   POWER DENSITY NORMAL TO BEAM
C     RADPOW(7:8,..)   NUMBER OF PHOTONS
C     RADPOW(9:10,..)  NUMBER OF PHOTONS NORMAL TO BEAM
C     RADPOW(11:12,..) 1D POWER DISTRIBUTION
C
C     RADPOW(13,..)    Z-COORDINATE ON ABSORBER
C     RADPOW(14,..)   DUMMY
C     RADPOW(15,..)    POWER DENSITY ON ABSORBER
C     RADPOW(16,..)   DUMMY
C     RADPOW(17,..)    NUMBER OF PHOTONS ON ABSORBER
C     RADPOW(18,..)   DUMMY
C     RADPOW(19,..)   1D POWER DENSITY ON ABSORBER
C LOGBUCH SEITE 23

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,track.
      include 'track.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER IPOL,NPOL,IPOI,INSIDE,IWALL,ID,IBIN,NBIN
      INTEGER IMODE,NN,IEND,ISTART,ITEST,IWARN1,IWARN2,IWARN3,IWARN4
      INTEGER ICYCLE
      DOUBLE PRECISION BYA,BY,Z,X,ZP,PX,PZ,DIS2,DXBIN,ZPSIGN,D2POW,SINPHI,COSPHI
      DOUBLE PRECISION D1POW,ECGAM,DNGAM,DZBIN,DISMIN,DIS1,DISMINA
      REAL*4 XI,XE,X1,X2,BYOLD
      real*8 hmaxm,hsumm
      REAL*4 XFILL,YFILL

C21.9.92     DOUBLE PRECISION CONS

      IF (POWBCUT.GE.0.0D0.AND.POWBCUT.LE.1.E-4) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** WARNING BEAMPOW ***'
           WRITE(LUNGFO,*)'POWBCUT.GE.0.0D0 .AND. POWBCUT.LE.1.E-4'
           WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
           WRITE(LUNGFO,*)
           WRITE(6,*)
           WRITE(6,*)'*** WARNING EAMPOW ***'
           WRITE(6,*)'POWBCUT.GE.0.0D0 .AND. POWBCUT.LE.1.E-4'
           WRITE(6,*)'CHECK RESULTS CAREFULLY'
           WRITE(6,*)
      ENDIF

      IF (IAMPLI.NE.0) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'*** WARNING IN BEAMPOW ***'
         WRITE(LUNGFO,*)'IAMPLI.NE.0, BUT IAMPLI NOT TAKEN INTO ACCOUNT'
         WRITE(LUNGFO,*)'FOR POWER CALCULATIONS WITH OPTION IPOWER'
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)'*** WARNING IN BEAMPOW ***'
         WRITE(6,*)'IAMPLI.NE.0, BUT IAMPLI NOT TAKEN INTO ACCOUNT'
         WRITE(6,*)'FOR POWER CALCULATIONS WITH OPTION IPOWER'
         WRITE(6,*)
      ENDIF !IAMPLI

      ALLOCATE(POWS2(NCO))
      ALLOCATE(RADPOW(19,2*NPOLMX,NCO))
      ALLOCATE(IPOLLIM(2,2*NPOLMX))

      POWCOR=1.   !3.5.93

      DISMIN=1.D30
      DISMINA=1.D30

      IF (XABSORB   .EQ.9999.) XABSORB=XWALLE
      IF (ZABSORB(1).EQ.9999.) ZABSORB(1)=WALL(1)
      IF (ZABSORB(2).EQ.9999.) ZABSORB(2)=WALL(2)

      IF (ZABSORB(1).LE.ZABSORB(2)) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** ERROR IN BEAMPOW ***'
           WRITE(LUNGFO,*)'Z-POSITION OF ABSORBER MEANINGLESS'
           WRITE(LUNGFO,*)'CHECK ZABSORB IN NAMELIST SPECTN'
           WRITE(LUNGFO,*)
           WRITE(6,*)
           WRITE(6,*)'*** ERROR IN BEAMPOW ***'
           WRITE(6,*)'Z-POSITION OF ABSORBER MEANINGLESS'
           WRITE(6,*)'CHECK ZABSORB IN NAMELIST SPECTN'
           WRITE(6,*)
           STOP
      ENDIF

      IF (XABSORB.LE.WTRA(1,1,NCO)) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** WARNING SR BEAMPOW ***'
           WRITE(LUNGFO,*)'ABSORBER BEFORE END OF TRAJECTORY'
           WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
           WRITE(LUNGFO,*)
           WRITE(6,*)
           WRITE(6,*)'*** WARNING SR BEAMPOW ***'
           WRITE(6,*)'ABSORBER BEFORE END OF TRAJECTORY'
           WRITE(6,*)'CHECK RESULTS CAREFULLY'
           WRITE(6,*)
      ENDIF

      IF (WALL(1).LE.WALL(2)) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** ERROR IN BEAMPOW ***'
           WRITE(LUNGFO,*)'BEAMLINE MEANINGLESS'
           WRITE(LUNGFO,*)'CHECK WALL IN NAMELIST SPECTN'
           WRITE(LUNGFO,*)
           WRITE(6,*)
           WRITE(6,*)'*** ERROR IN BEAMPOW ***'
           WRITE(6,*)'BEAMLINE MEANINGLESS'
           WRITE(6,*)'CHECK WALL IN NAMELIST SPECTN'
           WRITE(6,*)
           STOP
      ENDIF

      IF (XWALLE.LE.XWALLI) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** ERROR IN BEAMPOW ***'
           WRITE(LUNGFO,*)'BEAMLINE MEANINGLESS'
           WRITE(LUNGFO,*)'CHECK XWALLI,XWALLE IN NAMELIST SPECTN'
           WRITE(LUNGFO,*)
           WRITE(6,*)
           WRITE(6,*)'*** ERROR IN BEAMPOW ***'
           WRITE(6,*)'BEAMLINE MEANINGLESS'
           WRITE(6,*)'CHECK XWALLI,XWALLE IN NAMELIST SPECTN'
           WRITE(6,*)
           STOP
      ENDIF

      IF (WALL(1).LT.ZMX.OR.WALL(2).GT.ZMN) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** WARNING SR BEAMPOW ***'
           WRITE(LUNGFO,*)'WALLS OF BEAMLINE INCOMPATIBLE WITH TRAJECTORY'
           WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
           WRITE(LUNGFO,*)
           WRITE(6,*)
           WRITE(6,*)'*** WARNING SR BEAMPOW ***'
           WRITE(6,*)'WALLS OF BEAMLINE INCOMPATIBLE WITH TRAJECTORY'
           WRITE(6,*)
      ENDIF

C--- FIND POLES OF WLS

      NPOL=0
      INSIDE=0

C21.9.92  CONS=8.85D-5*CLIGHT1/(4.*PI1*EMASSG1)

      POWCOR=1.0  !C21.9.92

      DO IPOI=1,NCO

      BY=WTRA(2,3,IPOI)
      BYA=DABS(BY)
      ZP=WTRA(3,2,IPOI)/WTRA(1,2,IPOI)

      IF (INSIDE.EQ.0) THEN

          IF (BYA.GT.POWBCUT.AND.DABS(ZP).GT.1.D-10) THEN
         INSIDE=1
         NPOL=NPOL+1
         IF (NPOL.GT.NDPOL) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** ERROR IN BEAMPOW ***'
           WRITE(LUNGFO,*)'DIMENSION EXCEEDED NDPOLP'
           WRITE(LUNGFO,*)'CHECK NPOLMX AND POWBCUT IN NAMELIST SPECTN'
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'POLE NUMBER, (SLOPE ZP,X,BY) AT START, (SLOPE ZP,X,BY) AT END OF POLES ALREADY DETECTED:'
           DO IPOL=1,NPOL-1
             WRITE(LUNGFO,*)'POLES:',IPOL
             WRITE(LUNGFO,*)
     &              SNGL(WTRA(3,2,IPOLLIM(1,IPOL))
     &              /WTRA(1,2,IPOLLIM(1,IPOL)))
     &              ,SNGL(WTRA(1,1,IPOLLIM(1,IPOL)))
     &              ,SNGL(WTRA(2,3,IPOLLIM(1,IPOL)))
             WRITE(LUNGFO,*)
     &              SNGL(WTRA(3,2,IPOLLIM(2,IPOL))
     &              /WTRA(1,2,IPOLLIM(2,IPOL)))
     &              ,SNGL(WTRA(1,1,IPOLLIM(2,IPOL)))
     &              ,SNGL(WTRA(2,3,IPOLLIM(2,IPOL)))
           ENDDO  !NPOL
           WRITE(LUNGFO,*)'POLE NUMBER, (SLOPE ZP,X,BY) AT START OF CURRENT POLE:'
             WRITE(LUNGFO,*)'POLES:',NPOL
             WRITE(LUNGFO,*)
     &              SNGL(WTRA(3,2,IPOI)
     &              /WTRA(1,2,IPOI))
     &              ,SNGL(WTRA(1,1,IPOI))
     &              ,SNGL(WTRA(2,3,IPOI))
           WRITE(6,*)
           WRITE(6,*)'*** ERROR IN BEAMPOW ***'
           WRITE(6,*)'DIMENSION EXCEEDED NDPOLP'
           WRITE(6,*)'CHECK NPOLMX AND POWBCUT IN NAMELIST SPECTN'
           WRITE(6,*)
           WRITE(6,*)'POLE NUMBER, (SLOPE ZP,X,BY) AT START, (SLOPE ZP,X,BY) AT END OF POLES ALREADY DETECTED:'
           DO IPOL=1,NPOL-1
             WRITE(6,*)'POLES:',IPOL
             WRITE(6,*)
     &              SNGL(WTRA(3,2,IPOLLIM(1,IPOL))
     &              /WTRA(1,2,IPOLLIM(1,IPOL)))
     &              ,SNGL(WTRA(1,1,IPOLLIM(1,IPOL)))
     &              ,SNGL(WTRA(2,3,IPOLLIM(1,IPOL)))
             WRITE(6,*)
     &              SNGL(WTRA(3,2,IPOLLIM(2,IPOL))
     &              /WTRA(1,2,IPOLLIM(2,IPOL)))
     &              ,SNGL(WTRA(1,1,IPOLLIM(2,IPOL)))
     &              ,SNGL(WTRA(2,3,IPOLLIM(2,IPOL)))
           ENDDO  !NPOL
           WRITE(6,*)'POLE NUMBER, (SLOPE ZP,X,BY) AT START OF CURRENT POLE:'
             WRITE(6,*)'POLES:',NPOL
             WRITE(6,*)
     &              SNGL(WTRA(3,2,IPOI)
     &              /WTRA(1,2,IPOI))
     &              ,SNGL(WTRA(1,1,IPOI))
     &              ,SNGL(WTRA(2,3,IPOI))
           STOP
         ENDIF !NPOL
         ZPSIGN=DSIGN(1.D0,ZP)
         IPOLLIM(1,NPOL)=IPOI
          ENDIF   !POWBCUT

      ELSE  !INSIDE

          IF (BYOLD*BY.LE.0..OR.BYA.LT.POWBCUT) THEN
         INSIDE=0
C20.10.92      IPOLLIM(2,NPOL)=IPOI
         IPOLLIM(2,NPOL)=IPOI-1
          ENDIF   !POWBCUT)

C--- TRAJECTORY REACHES MAXIMAL DISPLACEMENT, I.E. LIGHT HITS OPPOSITE WALL

          IF(DSIGN(1.D0,ZP).NE.ZPSIGN) THEN
C20.10.92      IPOLLIM(2,NPOL)=IPOI
         IPOLLIM(2,NPOL)=IPOI-1
            IF (BYA.GT.POWBCUT.AND.DABS(ZP).GT.1.D-10) THEN
         NPOL=NPOL+1
         IF (NPOL.GT.NDPOL) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** ERROR IN BEAMPOW ***'
           WRITE(LUNGFO,*)'DIMENSION EXCEEDED NDPOLP'
           WRITE(LUNGFO,*)'CHECK NPOLMX AND POWBCUT IN NAMELIST SPECTN'
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'POLE NUMBER, (SLOPE ZP,X,BY) AT START, (SLOPE ZP,X,BY) AT END OF POLES ALREADY DETECTED:'
           DO IPOL=1,NPOL-1
             WRITE(LUNGFO,*)'POLES:',IPOL
             WRITE(LUNGFO,*)
     &              SNGL(WTRA(3,2,IPOLLIM(1,IPOL))
     &              /WTRA(1,2,IPOLLIM(1,IPOL)))
     &              ,SNGL(WTRA(1,1,IPOLLIM(1,IPOL)))
     &              ,SNGL(WTRA(2,3,IPOLLIM(1,IPOL)))
             WRITE(LUNGFO,*)
     &              SNGL(WTRA(3,2,IPOLLIM(2,IPOL))
     &              /WTRA(1,2,IPOLLIM(2,IPOL)))
     &              ,SNGL(WTRA(1,1,IPOLLIM(2,IPOL)))
     &              ,SNGL(WTRA(2,3,IPOLLIM(2,IPOL)))
           ENDDO  !IPOL
           WRITE(LUNGFO,*)'POLE NUMBER, (SLOPE ZP,X,BY) AT START OF CURRENT POLE:'
             WRITE(LUNGFO,*)'POLES:',NPOL
             WRITE(LUNGFO,*)
     &              SNGL(WTRA(3,2,IPOI)
     &              /WTRA(1,2,IPOI))
     &              ,SNGL(WTRA(1,1,IPOI))
     &              ,SNGL(WTRA(2,3,IPOI))
           WRITE(6,*)
           WRITE(6,*)'*** ERROR IN BEAMPOW ***'
           WRITE(6,*)'DIMENSION EXCEEDED NDPOLP'
           WRITE(6,*)'CHECK NPOLMX AND POWBCUT IN NAMELIST SPECTN'
           WRITE(6,*)
           WRITE(6,*)'POL NUMBER, (SLOPE ZP,X,BY) AT START, (SLOPE ZP,X,BY) AT END OF POLES ALREADY DETECTED:'
           DO IPOL=1,NPOL-1
             WRITE(6,*)'POLES:',IPOL
             WRITE(6,*)
     &              SNGL(WTRA(3,2,IPOLLIM(1,IPOL))
     &              /WTRA(1,2,IPOLLIM(1,IPOL)))
     &              ,SNGL(WTRA(1,1,IPOLLIM(1,IPOL)))
     &              ,SNGL(WTRA(2,3,IPOLLIM(1,IPOL)))
             WRITE(6,*)
     &              SNGL(WTRA(3,2,IPOLLIM(2,IPOL))
     &              /WTRA(1,2,IPOLLIM(2,IPOL)))
     &              ,SNGL(WTRA(1,1,IPOLLIM(2,IPOL)))
     &              ,SNGL(WTRA(2,3,IPOLLIM(2,IPOL)))
           ENDDO  !IPOL
           WRITE(6,*)'POLE NUMBER, (SLOPE ZP,X,BY) AT START OF CURRENT POLE:'
             WRITE(6,*)'POLES:',NPOL
             WRITE(6,*)
     &              SNGL(WTRA(3,2,IPOI)
     &              /WTRA(1,2,IPOI))
     &              ,SNGL(WTRA(1,1,IPOI))
     &              ,SNGL(WTRA(2,3,IPOI))
           STOP
         ENDIF !NPOL
              IPOLLIM(1,NPOL)=IPOI
         ZPSIGN=DSIGN(1.D0,ZP)
            ELSE  !POWCUT
         INSIDE=0
            ENDIF !POWCUT
          ENDIF   !DSIGN(1.D0,ZP)

      ENDIF !INSIDE

      BYOLD=BY
      ENDDO !IPOI

      IF (INSIDE.EQ.1) IPOLLIM(2,NPOL)=NCO


C--- LOOP OVER ALL POLES

      DO IPOL=1,NPOL

C--- LOOP OVER POINTS OF REFERENCE ORBIT

      DO IPOI=IPOLLIM(1,IPOL),IPOLLIM(2,IPOL)

C--- FIND POINT P WHERE RADIATION HITS THE WALL (LOGBOOK S.21)
C    AND CALCULATE POWER DENSITY

      X=WTRA(1,1,IPOI)
      Z=WTRA(3,1,IPOI)
      ZP=WTRA(3,2,IPOI)/WTRA(1,2,IPOI)
      BYA=DABS(WTRA(2,3,IPOI))
      SINPHI=DSQRT(ZP**2/(1.+ZP**2))
      COSPHI=DSQRT(1./(1.+ZP**2))

C--- POWER ON BEAMLINE WALL AND NORMAL TO BEAM

      DO IWALL=1,2

          PZ=WALL(IWALL)
          PX=(PZ-Z)/ZP+X
          DIS2=(X-PX)**2+(Z-PZ)**2
          DIS1=DSQRT(DIS2)
            IF (PX.LT.X) THEN
         D2POW=0.0
          ELSE
C21.9.92         D2POW=CONS*DMYENERGY**4*DMYCUR*BYA/DIS2*POWCOR
              D2POW=10.84/2.*1.D6*DMYENERGY**4*DMYCUR*BYA/DIS2
            ENDIF
            IF (IWALL.EQ.1.AND.Z.GT.WALL(IWALL))D2POW=0.0
            IF (IWALL.EQ.2.AND.Z.LT.WALL(IWALL))D2POW=0.0
          IF (D2POW.GT.0.0.AND.DIS1.LT.DISMIN
     &         .AND.PX.GE.XWALLI.AND.PX.LE.XWALLE) DISMIN=DIS1
C21.9.92     D1POW=D2POW*2.*DIS1/DMYGAMMA !D2POW INTEGRATED OVER Y
          IF(D2POW.NE.0.) THEN
             D1POW=(CGAM1*CLIGHT1/(2.D0*PI1))
     &               *DMYENERGY**3*DMYCUR*BYA/DIS1 !D2POW INTEGRATED OVER Y
          ELSE
             D1POW=0.0
          ENDIF
          ECGAM=ecdipev1*DMYENERGY**2*BYA*ECHARGE1 !CRITICAL PHOTONENERGY (JOULE)
          IF(ECGAM.NE.0.0) THEN
C240593     DNGAM=3.25*D1POW/ECGAM/POWCOR
C           RATE OF PHOTONS PER UNIT LENGTH
         DNGAM=15.*DSQRT(3.D0)/8.*D1POW/ECGAM/POWCOR
          ELSE
         DNGAM=0.0
          ENDIF   !ECGAM
          RADPOW(IWALL,IPOL,IPOI)=PX
          RADPOW(IWALL+ 2,IPOL,IPOI)=D2POW*SINPHI
          RADPOW(IWALL+ 4,IPOL,IPOI)=D2POW
          RADPOW(IWALL+ 6,IPOL,IPOI)=DNGAM*SINPHI
          RADPOW(IWALL+ 8,IPOL,IPOI)=DNGAM
          RADPOW(IWALL+10,IPOL,IPOI)=D1POW*SINPHI/POWCOR

      ENDDO !IWALL

      ENDDO !IPOI
      ENDDO !NPOL

C--- ASCENDING ORDER IN ARRAY

      DO IWALL=1,2
      DO IPOL=1,NPOL

          ISTART=IPOLLIM(1,IPOL)
C20.10.92       IEND=IPOLLIM(2,IPOL)-1
          IEND=IPOLLIM(2,IPOL)
          ITEST=ISTART+(IEND-ISTART)/2

          IF (RADPOW(IWALL,IPOL,ITEST).GT.RADPOW(IWALL,IPOL,ITEST+1)) THEN

         DO IMODE=0,5
           DO IPOI=ISTART,IEND
             RADPOW(14,IPOL,IPOI)=
     &              RADPOW(IWALL+2*IMODE,IPOL,IEND-IPOI+ISTART)
           ENDDO  !IPOI
           DO IPOI=ISTART,IEND
             RADPOW(IWALL+2*IMODE,IPOL,IPOI)=
     &              RADPOW(14,IPOL,IPOI)
           ENDDO  !IPOI
         ENDDO !IMODE

          ENDIF

      ENDDO !IPOL
      ENDDO !IWALL

C--- CHECK SPACING

      DO IWALL=1,2
      DO IPOL=1,NPOL
      DO IPOI=IPOLLIM(1,IPOL),IPOLLIM(2,IPOL)-1
          X1=RADPOW(IWALL,IPOL,IPOI)
          X2=RADPOW(IWALL,IPOL,IPOI+1)
          IF(X1.GE.X2) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'*** WARNING SR BEAMPOW ***'
         WRITE(LUNGFO,*)'BAD SPACING OF POINTS ON BEAMLINE WALL OCCURED'
           WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
         WRITE(LUNGFO,*)'TRY OTHER VALUES OF MYINUM OR OTHER BEAMLIME'
         WRITE(LUNGFO,*)'OR OTHER FIELD CONFIGURATION OR ...'
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)'*** WARNING SR BEAMPOW ***'
         WRITE(6,*)'BAD SPACING OF POINTS ON BEAMLINE WALL OCCURED'
           WRITE(6,*)'CHECK RESULTS CAREFULLY'
         WRITE(6,*)'TRY OTHER VALUES OF MYINUM OR OTHER BEAMLIME'
         WRITE(6,*)'OR OTHER FIELD CONFIGURATION OR ...'
         WRITE(6,*)
C20.10.92      STOP
          ENDIF
      ENDDO   !IPOI
      ENDDO !IPOL
      ENDDO !IWALL

C--- INTERPOLATE POWER DENSITY DISTRIBUTION BY SPLINES AND SUM UP
C    CONTRIBUTIONS OF ALL POLES, STORE INFORMATION IN HISTOGRAMS

      NBIN=NPWALL
      DXBIN=(XWALLE-XWALLI)/(NBIN-1)
      XI=XWALLI-DXBIN/2.
      XE=XI+NBIN*DXBIN

      ID=IDPOWER+1
      call hbook1m(ID,'2D POW. DENS. ON 1. WALL',
     &              NBIN,XI,XE,VMX)
      ID=IDPOWER+2
      call hbook1m(ID,'2D POW. DENS. ON 2. WALL',
     &              NBIN,XI,XE,VMX)

      ID=IDPOWER+1001
      call hbook1m(ID,'2D POW. DENS. ON 1. WALL, NORMAL',
     &              NBIN,XI,XE,VMX)
      ID=IDPOWER+1002
      call hbook1m(ID,
     &       '2D POW. DENS. ON 2. WALL, NORMAL',
     &        NBIN,XI,XE,VMX)

      ID=IDPOWER+2001
      call hbook1m(ID,'PHOTON RATE ON 1. WALL',
     &              NBIN,XI,XE,VMX)
      ID=IDPOWER+2002
      call hbook1m(ID,'PHOTON RATE ON 2. WALL',
     &              NBIN,XI,XE,VMX)

      ID=IDPOWER+3001
      call hbook1m(ID,'PHOTON RATE ON 1. WALL, NORMAL',
     &              NBIN,XI,XE,VMX)
      ID=IDPOWER+3002
      call hbook1m(ID,
     &       'PHOTON RATE ON 2. WALL, NORMAL',
     &        NBIN,XI,XE,VMX)

      ID=IDPOWER+4001
      call hbook1m(ID,'1D POW. DENS.,1. WALL',
     &              NBIN,XI,XE,VMX)
      ID=IDPOWER+4002
      call hbook1m(ID,'1D POW. DENS.,2. WALL',
     &              NBIN,XI,XE,VMX)


      DO IWALL=1,2
      DO IMODE=0,4

      ID=IDPOWER+1000*IMODE+IWALL

      DO IPOL=1,NPOL
C20.10.92       NN=IPOLLIM(2,IPOL)-IPOLLIM(1,IPOL)
          NN=IPOLLIM(2,IPOL)-IPOLLIM(1,IPOL)+1
      DO IBIN=1,NBIN

          XFILL=XI-DXBIN/2.+IBIN*DXBIN
          CALL POWINT(XFILL,YFILL,IWALL,IMODE,IPOL,
     &           IPOLLIM(1,IPOL),NN) !INTERPOLATION OF POWERDENSITY
          IF (YFILL.LT.0.AND.IWARN1.NE.1) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'*** WARNING SR BEAMPOW ***'
         WRITE(LUNGFO,*)'PROBLEMS WITH SPLINE-INTERPOLATION'
         WRITE(LUNGFO,*)'CHANGE SPACING'
         WRITE(LUNGFO,*)'NEGATIVE INTERPOLATION RESULT SET TO ZERO'
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)'*** WARNING SR BEAMPOW ***'
         WRITE(6,*)'PROBLEMS WITH SPLINE-INTERPOLATION'
         WRITE(6,*)'CHANGE SPACING'
         WRITE(6,*)'NEGATIVE INTERPOLATION RESULT SET TO ZERO'
         WRITE(6,*)
         IWARN1=1
          ENDIF !YFILL
          IF (YFILL.LT.0.) THEN
         YFILL=0.
          ENDIF !YFILL

          CALL hfillm(ID,XFILL,0.,dble(YFILL))

      ENDDO !IBIN
      ENDDO !IPOL
      ENDDO !IMODE
      ENDDO !IWALL

C--- TOTAL PHOTO DESORPTION

C DONT SUM UP NORMAL FLUX. RESULT WOULD BE WRONG SINCE
C INTEGRATION IS PERFORM ALONG THE BEAMLINE WHILE FLUX IS
C CALCULATED NORMAL TO BEAM

      DO IWALL=1,2
          TOTGAM(IWALL)  =HSUMM(IDPOWER+2000+IWALL)*DXBIN
          TOTGAM(IWALL+2)=HSUMM(IDPOWER+4000+IWALL)*DXBIN

          TOTMAX(IWALL)  =hmaxm(IDPOWER+0000+IWALL)
          TOTMAX(IWALL+2)=hmaxm(IDPOWER+3000+IWALL)
          TOTMAX(IWALL+4)=hmaxm(IDPOWER+4000+IWALL)
      ENDDO !IWALL

          TOTMAX(11)=DISMIN

C****************************************************************
C     NOW EVERYTHING FOR NORMAL ABSORBER
C****************************************************************

C--- FIND POLES OF WLS

      NPOL=0
      INSIDE=0

C21.9.92  CONS=8.85D-5*CLIGHT1/(4.*PI1*EMASSG1)

      POWCOR=1.0  !C21.9.92

      DO IPOI=1,NCO

      IF (WTRA(1,1,IPOI).GE.XIANF.AND.WTRA(1,1,IPOI).LE.XIEND) THEN
           BY=WTRA(2,3,IPOI)
      ELSE
           BY=0.0
      ENDIF


      BYA=DABS(BY)
      ZP=WTRA(3,2,IPOI)/WTRA(1,2,IPOI)

      IF (INSIDE.EQ.0) THEN

          IF (BYA.GT.POWBCUT) THEN
         INSIDE=1
         NPOL=NPOL+1
         IF (NPOL.GT.NDPOL) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** ERROR IN BEAMPOW ***'
           WRITE(LUNGFO,*)'DIMENSION EXCEEDED NDPOLP'
           WRITE(LUNGFO,*)'CHECK NPOLMX AND POWBCUT IN NAMELIST SPECTN'
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'POLE NUMBER, (SLOPE ZP,X,BY) AT START, (SLOPE ZP,X,BY) AT END OF POLES ALREADY DETECTED:'
           DO IPOL=1,NPOL-1
             WRITE(LUNGFO,*)'POLES:',IPOL
             WRITE(LUNGFO,*)
     &              SNGL(WTRA(3,2,IPOLLIM(1,IPOL))
     &              /WTRA(1,2,IPOLLIM(1,IPOL)))
     &              ,SNGL(WTRA(1,1,IPOLLIM(1,IPOL)))
     &              ,SNGL(WTRA(2,3,IPOLLIM(1,IPOL)))
             WRITE(LUNGFO,*)
     &              SNGL(WTRA(3,2,IPOLLIM(2,IPOL))
     &              /WTRA(1,2,IPOLLIM(2,IPOL)))
     &              ,SNGL(WTRA(1,1,IPOLLIM(2,IPOL)))
     &              ,SNGL(WTRA(2,3,IPOLLIM(2,IPOL)))
           ENDDO  !NPOL
           WRITE(LUNGFO,*)'POLE NUMBER, (SLOPE ZP,X,BY) AT START OF CURRENT POLE:'
             WRITE(LUNGFO,*)'POLES:',NPOL
             WRITE(LUNGFO,*)
     &              SNGL(WTRA(3,2,IPOI)
     &              /WTRA(1,2,IPOI))
     &              ,SNGL(WTRA(1,1,IPOI))
     &              ,SNGL(WTRA(2,3,IPOI))
           WRITE(6,*)
           WRITE(6,*)'*** ERROR IN BEAMPOW ***'
           WRITE(6,*)'DIMENSION EXCEEDED NDPOLP'
           WRITE(6,*)'CHECK NPOLMX AND POWBCUT IN NAMELIST SPECTN'
           WRITE(6,*)
           WRITE(6,*)'POLE NUMBER, (SLOPE ZP,X,BY) AT START, (SLOPE ZP,X,BY) AT END OF POLES ALREADY DETECTED:'
           DO IPOL=1,NPOL-1
             WRITE(6,*)'POLES:',IPOL
             WRITE(6,*)
     &              SNGL(WTRA(3,2,IPOLLIM(1,IPOL))
     &              /WTRA(1,2,IPOLLIM(1,IPOL)))
     &              ,SNGL(WTRA(1,1,IPOLLIM(1,IPOL)))
     &              ,SNGL(WTRA(2,3,IPOLLIM(1,IPOL)))
             WRITE(6,*)
     &              SNGL(WTRA(3,2,IPOLLIM(2,IPOL))
     &              /WTRA(1,2,IPOLLIM(2,IPOL)))
     &              ,SNGL(WTRA(1,1,IPOLLIM(2,IPOL)))
     &              ,SNGL(WTRA(2,3,IPOLLIM(2,IPOL)))
           ENDDO  !NPOL
           WRITE(6,*)'POLE NUMBER, (SLOPE ZP,X,BY) AT START OF CURRENT POLE:'
             WRITE(6,*)'POLES:',NPOL
             WRITE(6,*)
     &              SNGL(WTRA(3,2,IPOI)
     &              /WTRA(1,2,IPOI))
     &              ,SNGL(WTRA(1,1,IPOI))
     &              ,SNGL(WTRA(2,3,IPOI))
           STOP
         ENDIF !NPOL
         IPOLLIM(1,NPOL)=IPOI
          ENDIF   !POWBCUT

      ELSE  !INSIDE

          IF (BYOLD*BY.LE.0..OR.BYA.LT.POWBCUT) THEN
         INSIDE=0
C20.10.92      IPOLLIM(2,NPOL)=IPOI
         IPOLLIM(2,NPOL)=IPOI-1
          ENDIF   !POWBCUT

      ENDIF !INSIDE

      BYOLD=BY
      ENDDO !IPOI

      IF (INSIDE.EQ.1) IPOLLIM(2,NPOL)=NCO

C--- LOOP OVER ALL POLES

      DO IPOL=1,NPOL

C--- LOOP OVER POINTS OF REFERENCE ORBIT

      DO IPOI=IPOLLIM(1,IPOL),IPOLLIM(2,IPOL)

C--- FIND POINT P WHERE RADIATION HITS THE WALL (LOGBOOK S.21)
C    AND CALCULATE POWER DENSITY

      X=WTRA(1,1,IPOI)
      Z=WTRA(3,1,IPOI)
      ZP=WTRA(3,2,IPOI)/WTRA(1,2,IPOI)
      BYA=DABS(WTRA(2,3,IPOI))
      SINPHI=DSQRT(ZP**2/(1.+ZP**2))
      COSPHI=DSQRT(1./(1.+ZP**2))

C--- POWER

          PZ=Z+ZP*(XABSORB-X)
          PX=XABSORB
          DIS2=(X-PX)**2+(Z-PZ)**2
          DIS1=DSQRT(DIS2)
            IF (PX.LT.X) THEN
         D2POW=0.0
          ELSE
C21.9.92         D2POW=CONS*DMYENERGY**4*DMYCUR*BYA/DIS2*POWCOR
              D2POW=10.84/2.*1.D6*DMYENERGY**4*DMYCUR*BYA/DIS2
            ENDIF
          IF (D2POW.GT.0.0.AND.DIS1.LT.DISMINA) DISMINA=DIS1
C21.9.92     D1POW=D2POW*2.*DIS1/DMYGAMMA !D2POW INTEGRATED OVER Y
          IF(D2POW.NE.0.) THEN
            D1POW=(CGAM1*CLIGHT1/(2.D0*PI1))
     &              *DMYENERGY**3*DMYCUR*BYA/DIS1  !D2POW INTEGRATED OVER Y
          ELSE
            D1POW=0.0
          ENDIF
          ECGAM=ecdipev1*DMYENERGY**2*BYA*ECHARGE1 !CRITICAL PHOTONENERGY (JOULE)
          IF(ECGAM.NE.0.0) THEN
C240593     DNGAM=3.25*D1POW/ECGAM/POWCOR
C           RATE OF PHOTONS PER UNIT LENGTH
         DNGAM=15.*DSQRT(3.D0)/8.*D1POW/ECGAM/POWCOR
          ELSE
         DNGAM=0.0
          ENDIF   !ECGAM
          RADPOW(13,IPOL,IPOI)=PZ
          RADPOW(15,IPOL,IPOI)=D2POW*COSPHI
          RADPOW(17,IPOL,IPOI)=DNGAM*COSPHI
          RADPOW(19,IPOL,IPOI)=D1POW*COSPHI/POWCOR

      ENDDO !IPOI
      ENDDO !NPOL

C- ASCENDING ORDER

      DO IPOL=1,NPOL

          ISTART=IPOLLIM(1,IPOL)
          IEND=IPOLLIM(2,IPOL)

          IF (RADPOW(13,IPOL,ISTART).GT.RADPOW(13,IPOL,IEND)) THEN

         DO IMODE=0,3
           DO IPOI=ISTART,IEND
             RADPOW(14,IPOL,IPOI)=
     &              RADPOW(13+2*IMODE,IPOL,IEND-IPOI+ISTART)
           ENDDO  !IPOI
           DO IPOI=ISTART,IEND
             RADPOW(13+2*IMODE,IPOL,IPOI)=
     &              RADPOW(14,IPOL,IPOI)
           ENDDO  !IPOI
         ENDDO !IMODE

          ENDIF

      ENDDO !IPOL

C--- CHECK SPACING

      DO IPOL=1,NPOL
      DO IPOI=IPOLLIM(1,IPOL),IPOLLIM(2,IPOL)-1
          X1=RADPOW(13,IPOL,IPOI)
          X2=RADPOW(13,IPOL,IPOI+1)
          IF(X1.GE.X2) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'*** WARNING SR BEAMPOW ***'
         WRITE(LUNGFO,*)'BAD SPACING OF POINTS ON ABSORBER OCCURED'
           WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
         WRITE(LUNGFO,*)'TRY OTHER VALUES OF MYINUM OR ABSORBER'
         WRITE(LUNGFO,*)'OR OTHER FIELD CONFIGURATION OR ...'
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)'*** WARNING SR BEAMPOW ***'
         WRITE(6,*)'BAD SPACING OF POINTS ON ABSORBER OCCURED'
           WRITE(6,*)'CHECK RESULTS CAREFULLY'
         WRITE(6,*)'TRY OTHER VALUES OF MYINUM OR OTHER ABSORBER'
         WRITE(6,*)'OR OTHER FIELD CONFIGURATION OR ...'
         WRITE(6,*)
C20.10.92      STOP
          ENDIF
      ENDDO   !IPOI
      ENDDO !IPOL

C--- INTERPOLATE POWER DENSITY DISTRIBUTION BY SPLINES AND SUM UP
C    CONTRIBUTIONS OF ALL POLES, STORE INFORMATION IN HISTOGRAMS

      DZBIN=(ZABSORB(1)-ZABSORB(2))/(NBIN-1)
      XI=ZABSORB(2)-DZBIN/2.
      XE=XI+NBIN*DZBIN

      ID=IDPOWER+5001
      call hbook1m(ID,'2D POW. DENS. ON ABSORBER',
     &              NBIN,XI,XE,VMX)

      ID=IDPOWER+5002
      call hbook1m(ID,'PHOTON RATE ON ABSORBER',
     &              NBIN,XI,XE,VMX)

      ID=IDPOWER+5003
      call hbook1m(ID,'1D POW. DENS.,ABSORBER',
     &              NBIN,XI,XE,VMX)

      DO IPOL=1,NPOL
          NN=IPOLLIM(2,IPOL)-IPOLLIM(1,IPOL)+1
      DO IBIN=1,NBIN

          ID=IDPOWER+5001
          XFILL=XI-DZBIN/2.+IBIN*DZBIN
          CALL POWINT(XFILL,YFILL,13,0,IPOL,
     &           IPOLLIM(1,IPOL),NN) !INTERPOLATION OF POWERDENSITY
          IF (YFILL.LT.0.AND.IWARN2.NE.1) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'*** WARNING SR BEAMPOW ***'
         WRITE(LUNGFO,*)'PROBLEMS WITH SPLINE-INTERPOLATION'
         WRITE(LUNGFO,*)'CHANGE SPACING'
         WRITE(LUNGFO,*)'NEGATIVE INTERPOLATION RESULT SET TO ZERO'
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)'*** WARNING SR BEAMPOW ***'
         WRITE(6,*)'PROBLEMS WITH SPLINE-INTERPOLATION'
         WRITE(6,*)'CHANGE SPACING'
         WRITE(6,*)'NEGATIVE INTERPOLATION RESULT SET TO ZERO'
         WRITE(6,*)
         IWARN2=1
          ENDIF !YFILL
          IF (YFILL.LT.0.) THEN
         YFILL=0.
          ENDIF !YFILL

          CALL hfillm(ID,XFILL,0.,dble(YFILL))

          ID=IDPOWER+5002
          CALL POWINT(XFILL,YFILL,13,1,IPOL,
     &           IPOLLIM(1,IPOL),NN)
          IF (YFILL.LT.0.AND.IWARN3.NE.1) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'*** WARNING SR BEAMPOW ***'
         WRITE(LUNGFO,*)'PROBLEMS WITH SPLINE-INTERPOLATION'
         WRITE(LUNGFO,*)'CHANGE SPACING'
         WRITE(LUNGFO,*)'NEGATIVE INTERPOLATION RESULT SET TO ZERO'
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)'*** WARNING SR BEAMPOW ***'
         WRITE(6,*)'PROBLEMS WITH SPLINE-INTERPOLATION'
         WRITE(6,*)'CHANGE SPACING'
         WRITE(6,*)'NEGATIVE INTERPOLATION RESULT SET TO ZERO'
         WRITE(6,*)
         IWARN3=1
          ENDIF !YFILL
          IF (YFILL.LT.0.) THEN
         YFILL=0.
          ENDIF !YFILL

          CALL hfillm(ID,XFILL,0.,dble(yfill))

          ID=IDPOWER+5003
          CALL POWINT(XFILL,YFILL,13,2,IPOL,
     &           IPOLLIM(1,IPOL),NN)
          IF (YFILL.LT.0.AND.IWARN4.NE.1) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'*** WARNING SR BEAMPOW ***'
         WRITE(LUNGFO,*)'PROBLEMS WITH SPLINE-INTERPOLATION'
         WRITE(LUNGFO,*)'CHANGE SPACING'
         WRITE(LUNGFO,*)'NEGATIVE INTERPOLATION RESULT SET TO ZERO'
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)'*** WARNING SR BEAMPOW ***'
         WRITE(6,*)'PROBLEMS WITH SPLINE-INTERPOLATION'
         WRITE(6,*)'CHANGE SPACING'
         WRITE(6,*)'NEGATIVE INTERPOLATION RESULT SET TO ZERO'
         WRITE(6,*)
         IWARN4=1
          ENDIF !YFILL
          IF (YFILL.LT.0.) THEN
         YFILL=0.
          ENDIF !YFILL

          CALL hfillm(ID,XFILL,0.,dble(yfill))

      ENDDO !IBIN
      ENDDO !IPOL

      TOTGAM(5)=HSUMM(IDPOWER+5002)*DZBIN
      TOTGAM(6)=HSUMM(IDPOWER+5003)*DZBIN

      TOTMAX(7)=hmaxm(IDPOWER+5001)
      TOTMAX(8)=hmaxm(IDPOWER+5002)
      TOTMAX(9)=hmaxm(IDPOWER+5003)

      TOTMAX(12)=DISMINA

      CALL MHROUT(IDPOWER+1,ICYCLE,' ')
      CALL MHROUT(IDPOWER+2,ICYCLE,' ')
      CALL MHROUT(IDPOWER+1001,ICYCLE,' ')
      CALL MHROUT(IDPOWER+1002,ICYCLE,' ')
      CALL MHROUT(IDPOWER+2001,ICYCLE,' ')
      CALL MHROUT(IDPOWER+2002,ICYCLE,' ')
      CALL MHROUT(IDPOWER+3001,ICYCLE,' ')
      CALL MHROUT(IDPOWER+3002,ICYCLE,' ')
      CALL MHROUT(IDPOWER+4001,ICYCLE,' ')
      CALL MHROUT(IDPOWER+4002,ICYCLE,' ')
      CALL MHROUT(IDPOWER+5001,ICYCLE,' ')
      CALL MHROUT(IDPOWER+5002,ICYCLE,' ')
      CALL MHROUT(IDPOWER+5003,ICYCLE,' ')

      CALL hdeletm(IDPOWER+1)
      CALL hdeletm(IDPOWER+2)
      CALL hdeletm(IDPOWER+1001)
      CALL hdeletm(IDPOWER+1002)
      CALL hdeletm(IDPOWER+2001)
      CALL hdeletm(IDPOWER+2002)
      CALL hdeletm(IDPOWER+3001)
      CALL hdeletm(IDPOWER+3002)
      CALL hdeletm(IDPOWER+4001)
      CALL hdeletm(IDPOWER+4002)
      CALL hdeletm(IDPOWER+5001)
      CALL hdeletm(IDPOWER+5002)
      CALL hdeletm(IDPOWER+5003)

      DEALLOCATE(POWS2)

      RETURN
      END
