*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.57/00 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  15.55.02  by  Michael Scheer
*CMZ :  2.33/09 10/05/2001  18.02.55  by  Michael Scheer
*CMZ :  2.16/08 24/10/2000  11.19.57  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.43.04  by  Michael Scheer
*CMZ : 00.02/04 25/02/97  17.37.15  by  Michael Scheer
*CMZ : 00.01/09 01/09/95  13.03.01  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.45.26  by  Michael Scheer
*CMZ : 00.00/06 29/04/94  21.43.32  by  Michael Scheer
*CMZ : 00.00/05 29/04/94  20.18.16  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.15  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE CIRCPOW(ISOUR,NSOURCE)
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

*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.

C---  SUBROUTINE TO INTEGRATE SPECPOW OVER CIRCULAR PINHOLE
C     NZ,NY ARE ARRAY DIMENSIONS
C       INTEGRATION IS DONE FOR MZ AND MY POINTS

      IMPLICIT NONE
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEND.

      INTEGER NZ,NY,IY,IZ,MZ,MY,KLM,KLP,KHM,KHP,IIZ,IIY,ISOUR,IWCIRCO
      INTEGER IWSOUR,IWSOUR1,IWSOUR2,IWSOUR3,IWSOUR4,IDUM,NSOURCE,ICAL

      DOUBLE PRECISION SUM,SUMP,DZ,FL,FH,SH,SL,A,B,DY
     &  ,SPECK(NDOBSVP)
     &  ,SZ(NDOBSVZP),SY(NDOBSVYP),SUMY(NDOBSVYP),SUMYP(NDOBSVYP)
     &  ,FZ(NDOBSVZP),OBSVZF(NDOBSVZP),ZH,ZL
     &  ,DSUM,DSUML,DSUMH

      DATA ICAL/0/

      DO IDUM=1,NOBSV
        SPECK(IDUM)=SPECPOW(ISOUR+NSOURCE*(IDUM-1))
      ENDDO

      MZ=MOBSVZ
      NZ=NOBSVZ
      MY=MOBSVY
      NY=NOBSVY

      IF (IF1DIM.NE.0) THEN
        CALL CIRCPIN1D(NY,MY,SPECK,SUM,SUMP)
      SPECPOWVH(ISOUR)=SUM
        RETURN
      ENDIF

      IWSOUR=0
      IWCIRCO=IW_CIRC

      IF (MZ.LT.2) THEN
          WRITE(6,*) '*** ERROR IN CIRCPOW: MZ TOO LOW ***'
          STOP
      ENDIF
      IF (MY.LT.2) THEN
          WRITE(6,*) '*** ERROR IN CIRCPOW: MY TOO LOW ***'
          STOP
      ENDIF

C--- LOOP OF VERTICAL DIRECTION

      IIY=0
      DO IY=(NY-MY)/2+1,(NY-MY)/2+MY

      IIY=IIY+1
      SUMY(IIY)=0.D0

      IF (DABS(OBSVY(IY)-PINCEN(2)).LT.PINR-OBSVDY/4.D0) THEN

C--  SPLINE COEFFICENTS

      IIZ=0
      DO IZ=(NZ-MZ)/2+1,(NZ-MZ)+MZ
          IIZ=IIZ+1
          FZ(IIZ)=SPECK((IY-1)*NZ+IZ)
          OBSVZF(IIZ)=OBSVZ(IZ)
      ENDDO !IZ

      CALL FSPLINDX(OBSVDZ,FZ,IIZ,0.D0,0.D0,SZ)

C--  INTEGRATION LIMITS

      ZH=DSQRT(-(OBSVY(IY)-PINCEN(2))**2+PINR**2)
      ZL=PINCEN(3)-ZH
      ZH=PINCEN(3)+ZH

      CALL SPLINZY(IIZ,ZL,FL,OBSVZF,FZ,SZ,KLM)
      CALL SPLINZY(IIZ,ZH,FH,OBSVZF,FZ,SZ,KHM)

      KLP=KLM+1
      KHP=KHM+1

C--  CONTRIBUTION FORM PARTIAL INTERVALLS

      DZ=OBSVDZ
      DY=OBSVDY

      A=0.D0
      B=1.D0

      SH=
     &           -FZ(KLM)*A**2/2.*DZ+FZ(KLP)*B**2/2.*DZ
     &           +DZ**2/6.
     &           *((-DZ*A**4/4.+DZ*A**2/2.)*SZ(KLM)
     &           +(DZ*B**4/4.-DZ*B**2/2.)*SZ(KLP))

      A=(OBSVZF(KLP)-ZL)/OBSVDZ
      B=(ZL-OBSVZF(KLM))/OBSVDZ

      SL=
     &           -FZ(KLM)*A**2/2.*DZ+FZ(KLP)*B**2/2.*DZ
     &           +DZ**2/6.
     &           *((-DZ*A**4/4.+DZ*A**2/2.)*SZ(KLM)
     &           +(DZ*B**4/4.-DZ*B**2/2.)*SZ(KLP))

      DSUML=SH-SL

         IF(IWSOUR1.NE.ISOUR.AND.ISOUR.GT.0.AND.DSUML.LT.0.0) THEN
             IF ( DABS((SH-SL)/(SH+SL)) .GT. 5.D-15 ) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** WARNING IN CIRCPOW ***'
             WRITE(LUNGFO,*)
     &'SPLINE INTEGRATION FAILED FOR LEFT EDGE OF PINHOLE, RESULTS NOT RELIABLE'
             WRITE(LUNGFO,*)'SOURCE POINT:',ISOUR
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)
             WRITE(6,*)'*** WARNING IN CIRCPOW ***'
             WRITE(6,*)
     &'SPLINE INTEGRATION FAILED FOR LEFT EDGE OF PINHOLE, RESULTS NOT RELIABLE'
             WRITE(6,*)'SOURCE POINT:',ISOUR
             WRITE(6,*)
             WRITE(6,*)
             IWSOUR1=ISOUR
             IWSOUR=1
           IW_CIRC=1
             ENDIF    !SH-SL
         ENDIF !IWSOUR1

      A=1.D0
      B=0.D0

      SL=
     &           -FZ(KHM)*A**2/2.*DZ+FZ(KHP)*B**2/2.*DZ
     &           +DZ**2/6.
     &           *((-DZ*A**4/4.+DZ*A**2/2.)*SZ(KHM)
     &           +(DZ*B**4/4.-DZ*B**2/2.)*SZ(KHP))

      A=(OBSVZF(KHP)-ZH)/OBSVDZ
      B=(ZH-OBSVZF(KHM))/OBSVDZ

      SH=
     &           -FZ(KHM)*A**2/2.*DZ+FZ(KHP)*B**2/2.*DZ
     &           +DZ**2/6.
     &           *((-DZ*A**4/4.+DZ*A**2/2.)*SZ(KHM)
     &           +(DZ*B**4/4.-DZ*B**2/2.)*SZ(KHP))


      DSUMH=SH-SL

         IF(IWSOUR2.NE.ISOUR.AND.ISOUR.GT.0.AND.DSUMH.LT.0.0) THEN
             IF ( DABS((SH-SL)/(SH+SL)) .GT. 5.D-15 ) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** WARNING IN CIRCPOW ***'
             WRITE(LUNGFO,*)
     &          'SPLINE INTEGRATION FAILED FOR RIGHT EDGE OF PINHOLE, RESULTS NOT RELIABLE'
             WRITE(LUNGFO,*)'SOURCE POINT:',ISOUR
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)
             WRITE(6,*)'*** WARNING IN CIRCPOW ***'
             WRITE(6,*)
     &          'SPLINE INTEGRATION FAILED FOR RIGHT EDGE OF PINHOLE, RESULTS NOT RELIABLE'
             WRITE(6,*)'SOURCE POINT:',ISOUR
             WRITE(6,*)
             WRITE(6,*)
             IWSOUR2=ISOUR
             IWSOUR=1
           IW_CIRC=1
             ENDIF
         ENDIF !IWSOUR2

C---  NOW ALL FULL INTERVALLS

      DO IZ=KLP,KHM-1
         DSUM=
     &    OBSVDZ*0.5D0
     &    *(FZ(IZ)+FZ(IZ+1))
     &    -OBSVDZ**3/24.D0
     &    *(SZ(IZ)+SZ(IZ+1))

         IF(IWSOUR3.NE.ISOUR.AND.ISOUR.GT.0.AND.DSUM.LT.0.0) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** WARNING IN CIRCPOW ***'
             WRITE(LUNGFO,*)
     &          'SPLINE INTEGRATION FAILED INSIDE PINHOLE, RESULTS NOT RELIABLE'
             WRITE(LUNGFO,*)'SOURCE POINT:',ISOUR
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)
             WRITE(6,*)'*** WARNING IN CIRCPOW ***'
             WRITE(6,*)
     &              'SPLINE INTEGRATION FAILED INSIDE, RESULTS NOT RELIABLE'
             WRITE(6,*)'SOURCE POINT:',ISOUR
             WRITE(6,*)
             WRITE(6,*)
             IWSOUR3=ISOUR
             IWSOUR=1
           IW_CIRC=1
           DO IDUM=1,IIY
         SUMYP(IIY)=SUMY(IIY)
           ENDDO
         ENDIF !IWSOUR3

         SUMY(IIY)=SUMY(IIY)+DSUM
          IF (IW_CIRC.NE.0) SUMYP(IIY)=SUMYP(IIY)+ABS(DSUM)

      ENDDO !IZ

      SUMY(IIY)=SUMY(IIY)+DSUML+DSUMH
      IF (IW_CIRC.NE.0) SUMYP(IIY)=SUMYP(IIY)+ABS(DSUML)+ABS(DSUMH)

      ELSE !INSIDE PINHOLE

      SUMY(IIY)=0.0

      ENDIF !INSIDE PINHOLE
      ENDDO !IY

C--- INTEGRATION OVER VERTICAL DIRECTION

      CALL FSPLINDX(OBSVDY,SUMY,MY,0.D0,0.D0,SY)

      SUM=0.0
      DO IY=1,MY-1

        DSUM=
     &    OBSVDY*0.5D0
     &    *(SUMY(IY)+SUMY(IY+1))
     &    -OBSVDY**3/24.D0
     &    *(SY(IY)+SY(IY+1))

         IF(IWSOUR4.NE.ISOUR.AND.ISOUR.GT.0.AND.DSUM.LT.0.0) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** WARNING IN CIRCPOW ***'
             WRITE(LUNGFO,*)
     &          'SPLINE INTEGRATION FAILED FOR VERTICAL INTEGRATION, RESULTS NOT RELIABLE'
             WRITE(LUNGFO,*)'SOURCE POINT:',ISOUR
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)
             WRITE(6,*)'*** WARNING IN CIRCPOW ***'
             WRITE(6,*)
     &'SPLINE INTEGRATION FAILED FOR VERTICAL INTEGRATION, RESULTS NOT RELIABLE'
             WRITE(6,*)'SOURCE POINT:',ISOUR
             WRITE(6,*)
             WRITE(6,*)
             IWSOUR4=ISOUR
             IWSOUR=1
           IW_CIRC=1
         ENDIF !IWSOUR4

         SUM=SUM+DSUM

      ENDDO !IY

      IF (IW_CIRC.NE.0) THEN

      CALL FSPLINDX(OBSVDY,SUMYP,MY,0.D0,0.D0,SY)

      SUMP=0.0
      DO IY=1,MY-1

        DSUM=
     &    OBSVDY*0.5D0
     &    *(SUMYP(IY)+SUMYP(IY+1))
     &    -OBSVDY**3/24.D0
     &    *(SY(IY)+SY(IY+1))

         SUMP=SUMP+ABS(DSUM)

      ENDDO !IY

      ENDIF   !IW_CIRC

      SPECPOWVH(ISOUR)=SUM
      IW_CIRC=IWCIRCO

      IF (ICAL.EQ.0) THEN
      ICAL=1
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'*** WARNING IN CIRCPOW ***'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'FOR CIRCULAR PINHOLE ONLY POWER THROUGH PINHOLE'
      WRITE(LUNGFO,*)'IS CALCULATED. THE VERTICALLY INTEGRATED POWERDENSITY'
      WRITE(LUNGFO,*)'IS NOT, SINCE NOT MEANINGFUL'
      WRITE(LUNGFO,*)
      ENDIF

      RETURN
      END
