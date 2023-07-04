*CMZ :  4.00/13 02/09/2021  13.02.53  by  Michael Scheer
*CMZ :  4.00/11 10/05/2021  10.29.59  by  Michael Scheer
*CMZ :  3.04/00 18/01/2018  12.39.42  by  Michael Scheer
*CMZ :  3.01/00 16/07/2013  09.32.23  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.68/02 27/06/2012  16.34.34  by  Michael Scheer
*CMZ :  2.54/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.53/05 11/02/2005  09.55.20  by  Michael Scheer
*CMZ :  2.52/14 20/12/2004  17.10.56  by  Michael Scheer
*CMZ :  2.52/09 21/10/2004  15.47.48  by  Michael Scheer
*CMZ :  2.52/06 14/10/2004  09.16.20  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.34/07 04/09/2001  16.15.01  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.41.44  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  17.26.51  by  Michael Scheer
*CMZ : 00.01/11 11/09/96  17.24.24  by  Michael Scheer
*CMZ : 00.01/10 11/09/96  12.42.14  by  Michael Scheer
*CMZ : 00.01/07 16/03/95  14.21.07  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.45.58  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  18.05.04  by  Michael Scheer
*CMZ : 00.00/03 29/04/94  10.18.17  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.37  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BFOUR(XXIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,
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

C SUBROUTINE CALCULATES MAGNETIC FIELD AND VECTOR POTENTIAL FOR FOURIER
C EXPANSION OF HALBACH WIGGLER.
C FOURIER COEFFICIENTS ARE READ FROM DATA FILE OR CALCULATED
C ANALYTICALLY AND STORED IN COMMON-BLOCK FOR SIMPLE WAVELENGTH SHIFTER
C MODEL WITH ONE MAIN POLE.
C INTERNALLY HALBACHS CONVENTION IS USED FOR COORDINATE, INPUT AND
C OUTPUT CORRESPOND TO LAB.-SYSTEM.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,halbasy.
      include 'halbasy.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,fourier.
      include 'fourier.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      CHARACTER(60) CODEF

      INTEGER I,IK,ICAL,K,NKOEF,ICODEF,ifound

      DOUBLE PRECISION A0,A(MAXFOUR),XKFOUR(MAXFOUR),YKFOUR(MAXFOUR),ZKFOUR(MAXFOUR)
      REAL*4 AR(MAXFOUR),SUMCK2,SUMC
      DOUBLE PRECISION ZRFOUR,XL0FOUR,YL0FOUR,ZL0FOUR,XK0FOUR,YK0FOUR,ZK0FOUR,
     &  ZL0FOUR2,DUMZ,DUMN

      DOUBLE PRECISION XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT,
     &  DSNXKX,DCSXKX,DSHYKY,DCHYKY,DSNZKZ,DCSZKZ
     &  ,BXH,BYH,BZH,AXH,AYH,AZH,AN,AM,X,xxin


      DOUBLE PRECISION EXPOMY,DEXPOMY,EXPOMY1

      DOUBLE PRECISION DNULL
      COMPLEX*16 CDEXPOMX,CEXPOMZ,CDEXPOMZ
      DATA DNULL/0.D0/

      DATA ICAL/0/

      IF (KBEXTERN.NE.0) IRFILF=1

      xin=xxin+xshbfour

      IF (ICAL.NE.1) THEN

        IF (NFOUR.NE.-9999.AND.(NFOUR.LT.1.OR.NFOUR.GT.MAXFOUR)) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN BFOUR ***'
          WRITE(LUNGFO,*)'NFOUR.LT.1.OR.NFOUR.GT.MAXFOUR'
          WRITE(LUNGFO,*)'CHECK NFOUR (NAMELIST FOURIER)'
          WRITE(LUNGFO,*)'OR INCREASE MAXFOUR (FILE CMPARA.CMN)'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BFOUR ***'
          WRITE(6,*)'NFOUR.LT.1.OR.NFOUR.GT.MAXFOUR'
          WRITE(6,*)'CHECK NFOUR (NAMELIST FOURIER)'
          WRITE(6,*)'OR INCREASE MAXFOUR (FILE CMPARA.CMN)'
          WRITE(6,*)
          STOP
        ENDIF

        IF (IRFILF.NE.0) THEN

          OPEN(UNIT=LUNF,FILE=FILEF,STATUS='OLD',FORM='FORMATTED')

c          READ(LUNF,1000)ICODEF,CODEF
c1000      FORMAT(I10,'  ',1A60)
          read(lunf,'(a)')codef
          read(codef,*) icodef
          ifound=0
          do i=1,60
            if (
     &        codef(i:i).eq."0".or.
     &        codef(i:i).eq."1".or.
     &        codef(i:i).eq."2".or.
     &        codef(i:i).eq."3".or.
     &        codef(i:i).eq."4".or.
     &        codef(i:i).eq."5".or.
     &        codef(i:i).eq."6".or.
     &        codef(i:i).eq."7".or.
     &        codef(i:i).eq."8".or.
     &        codef(i:i).eq."9") ifound=1
            if (ifound.eq.1.and.codef(i:i).eq.' ') goto 1000
          enddo
1000      codef=codef(i+1:60)
          READ(LUNF,*)ZRFOUR
          READ(LUNF,*)NKOEF
          IF (NFOUR.EQ.-9999) NFOUR=NKOEF
          IF (NKOEF.GT.MAXFOUR.OR.NKOEF.LT.NFOUR) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
     &        '*** ERROR IN BFOUR ***'
            WRITE(LUNGFO,*)'NKOEF.GT.MAXFOUR.OR.NKOEF.LT.NFOUR'
            WRITE(LUNGFO,*)'CHECK PARAMETERS IN FILE CMPARA.CMN AND NAMELISTS FOURIER (AND MAYBE HALBASY)'
            WRITE(6,*)
            WRITE(6,*)
     &        '*** ERROR IN BFOUR ***'
            WRITE(6,*)'NKOEF.GT.MAXFOUR.OR.NKOEF.LT.NFOUR'
            WRITE(6,*)'CHECK PARAMETERS IN FILE CMPARA.CMN AND NAMELISTS FOURIER (AND MAYBE HALBASY)'
            STOP
          ENDIF

          READ(LUNF,*)IK,A0

          DO I=1,NKOEF-1
            READ(LUNF,*)IK,AR(IK-1)
            A(IK-1)=DBLE(AR(IK-1))
          END DO

          ZL0FOUR=ZRFOUR
          ZL0FOUR2=ZL0FOUR/2.
          ZK0FOUR=2.D0*PI1/ZL0FOUR
C021291          XL0FOUR=XLHALBASY  !VORERST SIEHE OBEN
          XL0FOUR=XLENFOUR      !VORERST SIEHE OBEN
          XK0FOUR=0.D0
          IF(XL0FOUR.NE.0.) XK0FOUR=2.D0*PI1/XL0FOUR
          YK0FOUR=DSQRT(ZK0FOUR**2+XK0FOUR**2)
          YL0FOUR=2.D0*PI1/YK0FOUR

          if (scbfour.eq.-9999.0d0) scbfour=1.0d0
          if (xshbfour.eq.-9999.0d0) xshbfour=0.0d0
          if (fouentr.eq.-9999.0d0) fouentr=-zl0four/2.0d0-xshbfour
          if (fouexit.eq.-9999.0d0) fouexit=zl0four/2.0d0-xshbfour

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     SR BFOUR:',NKOEF,' coefficients read from file:'
          WRITE(LUNGFO,*)'     ',FILEF
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     FOUENTR, FOUEXIT: ',fouentr,fouexit
          WRITE(LUNGFO,*)'     SCBFOUR, XSHBFOUR:',scbfour,xshbfour
          WRITE(LUNGFO,*)

          CLOSE(LUNF)

        ELSE !(IRFILF.NE.0)

          IF (AHWPOL.NE.1.) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN BFOUR ***'
            WRITE(LUNGFO,*)
     &        'FOR THIS FEATURE (IAHWPOL.NE.0) PARAMETER AHWPOL HAS TO BE 1. ***'
            WRITE(6,*)
            WRITE(6,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN BFOUR ***'
            WRITE(6,*)
     &        'FOR THIS FEATURE (IAHWPOL.NE.0) PARAMETER AHWPOL HAS TO BE 1. ***'
            WRITE(6,*)
            STOP
          ENDIF

          IF (NFOUR.EQ.-9999) NFOUR=NFOURWLS

          IF (FASYM.EQ.2.0D0) THEN
            ZL0FOUR=ZLHALBASY
          ELSE
            ZL0FOUR=ZLHALBASY/2.*(1.+FASYM)
          ENDIF

          ZL0FOUR2=ZL0FOUR/2.
          ZK0FOUR=2.D0*PI1/ZL0FOUR
C110996              XL0FOUR=XLHALBASY  !VORERST SIEHE OBEN
          XL0FOUR=XLENFOUR
          XK0FOUR=0.D0
          IF(XL0FOUR.NE.0.) XK0FOUR=2.D0*PI1/XL0FOUR
          YK0FOUR=DSQRT(ZK0FOUR**2+XK0FOUR**2)
          YL0FOUR=2.D0*PI1/YK0FOUR

          AN=FASYM

          A0=0.
          DO I=1,NFOUR-1

            AM=I

            DUMZ=
     &        -2.D0*
     &         (
     &          (2.D0*AM**2*AN**2-4.*AM**2-AN**2-2.D0*AN-1.D0)*
     &          DCOS((AM*PI1)/(AN+1.D0))-
     &          (2.D0*AM+AN+1.D0)*(2.D0*AM-AN-1.D0)*
     &          DCOS(AM*PI1)
     &         )*
     &         (AN+1.D0)

            DUMN=
     &        (AM*AN+AN+1.D0)*
     &        (AM*AN-AN-1.D0)*(2.D0*AM+AN+1.D0)*(2.D0*AM-AN-1.D0)*PI1

            IF (DUMN.NE.0.0D0) THEN
              A(I)=B0HALBASY*DUMZ/DUMN       !NORMIERUNG AUF PEAK-FELD
            ELSE
              A(I)=0.0D0
            ENDIF

          END DO

          WRITE(LUNGFO,*)
     &      '     SR BFOUR: Fourier expansion of wavelength shifter calculated'
          WRITE(LUNGFO,*)
     &      '               number of coefficients:',NFOUR
          WRITE(LUNGFO,*)

        ENDIF

        DO I=1,NFOUR-1
          ZKFOUR(I)=ZK0FOUR*I
          XKFOUR(I)=XK0FOUR !VORERST, SIEHE AUCH OBEN
          YKFOUR(I)=DSQRT(ZKFOUR(I)**2+XKFOUR(I)**2)
        END DO

C--- LONGITUDINAL FIELD HOMOGENITY, FIELD IS EXPANDED IN TAYLOR SERIE
C    TO SECOND ORDER. MAXIMUM AT Z=0 ASSUMED

CERROR 19MAY05
        SUMC=A0/2.0D0 !SUMC was not initialized (before 19may05)
        SUMCK2=0.0D0
        DO I=1,NFOUR-1
          SUMC=SUMC+A(I)
          SUMCK2=SUMCK2+A(I)*ZKFOUR(I)**2
        ENDDO
        IF (SUMCK2.NE.0.0) THEN
          XBHOMF=DSQRT(DABS(DBHOMF*2.D0/SUMCK2*SUMC))
        ELSE
          XBHOMF=9999.0D9
        ENDIF
C---

        IF (IRFILF.NE.0) THEN
          IF(IFOUR0.NE.0) A0=0.
          WRITE(LUNGFO,*)
     &      '     run number and comment of job that calculated Fourier coefficients:'
          WRITE(LUNGFO,*)ICODEF,'  ',CODEF
          WRITE(LUNGFO,*)'     flag IFOUR0:  ',IFOUR0
          WRITE(LUNGFO,*)'     period length ZL0FOUR:',SNGL(ZL0FOUR)
          WRITE(LUNGFO,*)
     &      '     lx, kx/kz:',SNGL(XL0FOUR),SNGL(XK0FOUR/ZK0FOUR)
          WRITE(LUNGFO,*)
        ENDIF

        IF(IPRNTF.NE.0) THEN
          WRITE(LUNGFO,*)'     Fourier coefficients and  kx/kz:'
          WRITE(LUNGFO,*)
          I=1
          WRITE(LUNGFO,*)'          0',SNGL(A0)
          DO I=1,NFOUR-1
            WRITE(LUNGFO,*)I,SNGL(A(I)),SNGL(XKFOUR(I)/ZKFOUR(I))
          ENDDO
        ENDIF !IPRNTF

        WRITE(LUNGFO,*)
     &    '     required homogeneity (DBHOMF) and corresponding distance from center:'
        WRITE(LUNGFO,*)'     ',SNGL(DBHOMF),SNGL(XBHOMF)
        WRITE(LUNGFO,*)

        DEVLEN=ZL0FOUR
        DEVLEN2=ZL0FOUR2

        ICAL=1

      ENDIF

C-------------------------------------------------------------------
C     IF (XIN.LT.-ZL0FOUR2.OR.XIN.GT.ZL0FOUR2) THEN  !VORSICHT WEGEN TRANPOLY
C     IF (DABS(XIN.LT.-ZL0FOUR2.OR.XIN.GT.ZL0FOUR2) THEN  !SIEHE LOGBUCH S. 218
C
      BXOUT=0.0
      BYOUT=0.0
      BZOUT=0.0
C
      AXOUT=0.
      AYOUT=0.
      AZOUT=0.
C
C        RETURN
C     ENDIF

      if (xin.lt.fouentr.or.xin.gt.fouexit) then
        return
      endif

      X=DMOD(XIN,ZL0FOUR) !2.12.91

      IF (X.GT.ZL0FOUR2) THEN
        X=X-ZL0FOUR
      ELSE IF (X.LT.-ZL0FOUR2) THEN
        X=X+ZL0FOUR
      ENDIF

      BXH=0.
      BYH=A0/2.D0
      BZH=0.

C IF CHANGED, CONSIDER FOLLOWING LOOP AND SR MYBFELD {



      AXH= A0/2.D0*  XIN !XIN IS HERE Z
      AYH=0.
      AZH=0.

C IF CHANGED, CONSIDER FOLLOWING LOOP AND SR MYBFELD }

      CDEXPOMX=CDEXP(DCMPLX(DNULL,XKFOUR(1)*(-ZIN)))
      DCSXKX=DREAL(CDEXPOMX)
      DSNXKX=DIMAG(CDEXPOMX)

      DEXPOMY=DEXP(YKFOUR(1)*YIN)
      EXPOMY=1.D0

      CDEXPOMZ=CDEXP(DCMPLX(DNULL,ZKFOUR(1)*    X ))
      CEXPOMZ=DCMPLX(1.D0,DNULL)

      DO K=1,NFOUR-1

        IF (XK0FOUR.NE.0.0D0) THEN
          EXPOMY=DEXP(YKFOUR(K)*YIN)
        ELSE
          EXPOMY=EXPOMY*DEXPOMY
        ENDIF

        EXPOMY1=1.D0/EXPOMY
        DCHYKY=(EXPOMY+EXPOMY1)*0.5D0
        DSHYKY=(EXPOMY-EXPOMY1)*0.5D0

        CEXPOMZ=CEXPOMZ*CDEXPOMZ
        DCSZKZ=DREAL(CEXPOMZ)
        DSNZKZ=DIMAG(CEXPOMZ)

        BXH=BXH-A(K)*XKFOUR(K)/YKFOUR(K)*DSNXKX*DSHYKY*DCSZKZ
        BYH=BYH+A(K)*                    DCSXKX*DCHYKY*DCSZKZ
        BZH=BZH-A(K)*ZKFOUR(K)/YKFOUR(K)*DCSXKX*DSHYKY*DSNZKZ



        AXH=AXH+A(K)/ZKFOUR(K)*DCSXKX*DCHYKY*DSNZKZ
        AZH=AZH+0.0

        AYH=AYH+A(K)/ZKFOUR(K)*XKFOUR(K)/YKFOUR(K)*DSNXKX*DSHYKY*DSNZKZ

      ENDDO

      BZOUT=-BXH*scbfour
      BYOUT= BYH*scbfour
      BXOUT= BZH*scbfour

      AZOUT=-AXH*scbfour
      AYOUT= AYH*scbfour
      AXOUT= AZH*scbfour

      IF (KBEXTERN.NE.0) IRFILF=0

      RETURN

      END
