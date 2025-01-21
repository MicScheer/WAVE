*CMZ :          21/01/2025  16.23.35  by  Michael Scheer
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
      SUBROUTINE BSEQFOUR(imag,XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)

      use magseqf90m

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
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      CHARACTER(60) CODEF

      INTEGER imag,I,IK,K,ICODEF,lunf

      DOUBLE PRECISION XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT,
     &  DSNXKX,DCSXKX,DSHYKY,DCHYKY,DSNZKZ,DCSZKZ
     &  ,BXH,BYH,BZH,AXH,AYH,AZH,AN,AM,X

      DOUBLE PRECISION EXPOMY,DEXPOMY,EXPOMY1

      DOUBLE PRECISION DNULL
      COMPLEX*16 CDEXPOMX,CEXPOMZ,CDEXPOMZ
      DATA DNULL/0.0D0/

      IF (seqmag(imag)%ical.NE.1) THEN

        OPEN(newUNIT=LUNF,FILE=seqmag(imag)%cfilefour,STATUS='OLD',FORM='FORMATTED')

        read(lunf,'(a)')seqmag(imag)%comm
        read(lunf,*) seqmag(imag)%zl0four
        READ(LUNF,*)seqmag(imag)%nfour
        allocate(
     &    seqmag(imag)%four(seqmag(imag)%nfour),
     &    seqmag(imag)%xkfour(seqmag(imag)%nfour),
     &    seqmag(imag)%ykfour(seqmag(imag)%nfour),
     &    seqmag(imag)%zkfour(seqmag(imag)%nfour)
     &    )

        READ(LUNF,*)IK,seqmag(imag)%A0

        DO I=1,seqmag(imag)%nfour-1
          READ(LUNF,*)IK,seqmag(imag)%four(IK-1)
        END DO

        seqmag(imag)%zl0four=seqmag(imag)%zl0four
        seqmag(imag)%zl0four2=seqmag(imag)%zl0four/2.
        seqmag(imag)%zk0four=2.D0*PI1/seqmag(imag)%zl0four
        seqmag(imag)%xl0four=seqmag(imag)%pmag(6)
        seqmag(imag)%xk0four=0.D0
        IF(seqmag(imag)%xl0four.NE.0.) seqmag(imag)%xk0four=2.D0*PI1/seqmag(imag)%xl0four
        seqmag(imag)%yk0four=DSQRT(seqmag(imag)%zk0four**2+seqmag(imag)%xk0four**2)
        seqmag(imag)%yl0four=2.D0*PI1/seqmag(imag)%yk0four

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     SR BSEQFOUR:',seqmag(imag)%nfour,' coefficients read from file:'
        WRITE(LUNGFO,*)'     ',seqmag(imag)%cfilefour
        WRITE(LUNGFO,*)

        CLOSE(LUNF)

        DO I=1,seqmag(imag)%nfour-1
          seqmag(imag)%zkfour(I)=seqmag(imag)%zk0four*I
          seqmag(imag)%xkfour(I)=seqmag(imag)%xk0four !VORERST, SIEHE AUCH OBEN
          seqmag(imag)%ykfour(I)=DSQRT(seqmag(imag)%zkfour(I)**2+seqmag(imag)%xkfour(I)**2)
        END DO

C--- LONGITUDINAL FIELD HOMOGENITY, FIELD IS EXPANDED IN TAYLOR SERIE
C    TO SECOND ORDER. MAXIMUM AT Z=0 ASSUMED

        WRITE(LUNGFO,*)
     &    '     Comment on file:'
        WRITE(LUNGFO,*)'     ',trim(seqmag(imag)%comm)
        WRITE(LUNGFO,*)'     period length seqmag(imag)%zl0four:',SNGL(seqmag(imag)%zl0four)
        WRITE(LUNGFO,*)
     &    '     lx, kx/kz:',SNGL(seqmag(imag)%xl0four),SNGL(seqmag(imag)%xk0four/seqmag(imag)%zk0four)
        WRITE(LUNGFO,*)

        DEVLEN=seqmag(imag)%zl0four
        DEVLEN2=seqmag(imag)%zl0four2

         seqmag(imag)%xcen= seqmag(imag)%pmag(1)
         seqmag(imag)%ycen= seqmag(imag)%pmag(2)
         seqmag(imag)%zcen= seqmag(imag)%pmag(3)

         seqmag(imag)%xmin= seqmag(imag)%xcen-devlen2
         seqmag(imag)%xmax= seqmag(imag)%xcen+devlen2

        seqmag(imag)%bmax=seqmag(imag)%A0/2.0D0
        DO I=1,seqmag(imag)%NFOUR-1
          seqmag(imag)%bmax=seqmag(imag)%bmax+seqmag(imag)%four(I)
        ENDDO

        seqmag(imag)%rho=seqmag(imag)%pmag(5)
        seqmag(imag)%scalebglob=-dbrho/seqmag(imag)%rho/seqmag(imag)%bmax

        seqmag(imag)%ical=1

      ENDIF

C-------------------------------------------------------------------
C     IF (XIN.LT.-seqmag(imag)%zl0four2.OR.XIN.GT.seqmag(imag)%zl0four2) THEN  !VORSICHT WEGEN TRANPOLY
C     IF (DABS(XIN.LT.-seqmag(imag)%zl0four2.OR.XIN.GT.seqmag(imag)%zl0four2) THEN  !SIEHE LOGBUCH S. 218
C
      BXOUT=0.0d0
      BYOUT=0.0d0
      BZOUT=0.0d0
C
      AXOUT=0.0d0
      AYOUT=0.0d0
      AZOUT=0.0d0
C
C        RETURN
C     ENDIF

      if (xin.lt.seqmag(imag)%xmin.or.xin.gt.seqmag(imag)%xmax) then
        return
      endif

      X=DMOD(XIN- seqmag(imag)%xcen,seqmag(imag)%zl0four) !2.12.91

      IF (X.GT.seqmag(imag)%zl0four2) THEN
        X=X-seqmag(imag)%zl0four
      ELSE IF (X.LT.-seqmag(imag)%zl0four2) THEN
        X=X+seqmag(imag)%zl0four
      ENDIF

      BXH=0.0d0
      BYH=seqmag(imag)%a0/2.D0
      BZH=0.0d0

C IF CHANGED, CONSIDER FOLLOWING LOOP AND SR MYBFELD {


      AXH= seqmag(imag)%a0/2.D0*  XIN !XIN IS HERE Z
      AYH=0.0d0
      AZH=0.0d0
C IF CHANGED, CONSIDER FOLLOWING LOOP AND SR MYBFELD }

      CDEXPOMX=CDEXP(DCMPLX(DNULL,seqmag(imag)%xkfour(1)*(-ZIN)))
      DCSXKX=DREAL(CDEXPOMX)
      DSNXKX=DIMAG(CDEXPOMX)

      DEXPOMY=DEXP(seqmag(imag)%ykfour(1)*YIN)
      EXPOMY=1.0D0

      CDEXPOMZ=CDEXP(DCMPLX(DNULL,seqmag(imag)%zkfour(1)*    X ))
      CEXPOMZ=DCMPLX(1.D0,DNULL)

      DO K=1,seqmag(imag)%nfour-1

        IF (seqmag(imag)%xk0four.NE.0.0D0) THEN
          EXPOMY=DEXP(seqmag(imag)%ykfour(K)*YIN)
        ELSE
          EXPOMY=EXPOMY*DEXPOMY
        ENDIF

        EXPOMY1=1.D0/EXPOMY
        DCHYKY=(EXPOMY+EXPOMY1)*0.5D0
        DSHYKY=(EXPOMY-EXPOMY1)*0.5D0

        CEXPOMZ=CEXPOMZ*CDEXPOMZ
        DCSZKZ=DREAL(CEXPOMZ)
        DSNZKZ=DIMAG(CEXPOMZ)

        BXH=BXH-seqmag(imag)%four(k)*seqmag(imag)%xkfour(K)/seqmag(imag)%ykfour(K)*DSNXKX*DSHYKY*DCSZKZ
        BYH=BYH+seqmag(imag)%four(k)*                    DCSXKX*DCHYKY*DCSZKZ
        BZH=BZH-seqmag(imag)%four(k)*seqmag(imag)%zkfour(K)/seqmag(imag)%ykfour(K)*DCSXKX*DSHYKY*DSNZKZ



        AXH=AXH+seqmag(imag)%four(k)/seqmag(imag)%zkfour(K)*DCSXKX*DCHYKY*DSNZKZ
        AZH=AZH+0.0

        AYH=AYH+seqmag(imag)%four(k)/seqmag(imag)%zkfour(K)*seqmag(imag)%xkfour(K)/seqmag(imag)%ykfour(K)*DSNXKX*DSHYKY*DSNZKZ

      ENDDO

      BZOUT=-BXH*seqmag(imag)%scalebglob
      BYOUT= BYH*seqmag(imag)%scalebglob
      BXOUT= BZH*seqmag(imag)%scalebglob

      AZOUT=-AXH*seqmag(imag)%scalebglob
      AYOUT= AYH*seqmag(imag)%scalebglob
      AXOUT= AZH*seqmag(imag)%scalebglob

      RETURN
      END
