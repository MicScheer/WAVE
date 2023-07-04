*CMZ :  2.70/05 02/01/2013  12.22.05  by  Michael Scheer
*CMZ :  2.65/02 28/09/2009  08.44.50  by  Michael Scheer
*CMZ :  2.63/04 11/06/2009  12.20.28  by  Michael Scheer
*CMZ :  2.48/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.47/23 17/02/2004  10.07.04  by  Michael Scheer
*CMZ :  2.47/18 26/11/2003  18.23.46  by  Michael Scheer
*CMZ :  2.47/08 08/05/2003  12.49.51  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.41/08 02/08/2002  19.53.38  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  14.39.20  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  17.16.55  by  Michael Scheer
*CMZ :  2.14/02 19/04/2000  17.02.45  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  11.45.40  by  Michael Scheer
*CMZ :  1.03/06 07/08/98  10.05.39  by  Michael Scheer
*CMZ : 00.02/04 12/02/97  11.00.58  by  Michael Scheer
*CMZ : 00.02/03 04/02/97  16.38.06  by  Michael Scheer
*-- Author :    Michael Scheer   22/01/97

      SUBROUTINE BPOLY2DH(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)
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

C *** ATTENTION: DEFINITION OF COEFFICIENTS DIFFERENT FROM R. WALKER'S

C     MAGNETIC FIELD ACCORDING TO POLYNOMIAL FIT
C     TRANSVERSALLY AND HARMONICS IN LONGITUDINALLY
C
C     INPUT/OUTPUT COORDINATE SYSTEM: X LONG., Y VERTICAL
C     INTERNAL COORDINATE SYSTEM: Z LONG., Y VERTICAL
C
C     UNIT: METER AND TESLA
C

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bpoly2dh.
      include 'bpoly2dh.cmn'
*KEND.

      INTEGER ICAL,LUNIN,IPAR,IHV,NN
      INTEGER N
      DOUBLE PRECISION N2,N4,N6,N8,N10

      DOUBLE PRECISION PI,XOFF
      DOUBLE PRECISION XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT
      DOUBLE PRECISION X,Y,Z,BXX,BYY,BZZ
      DOUBLE PRECISION XX,X2,X3,X4,X5,X6,X7,X8,X9,X10
      DOUBLE PRECISION YY,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10
      DOUBLE PRECISION XX2,XX3,XX4,XX5,XX6,XX7,XX8,XX9,XX10
      DOUBLE PRECISION YY2,YY3,YY4,YY5,YY6,YY7,YY8,YY9,YY10
      DOUBLE PRECISION ZL,ZK,ZK2,ZK4,ZK6,ZK8,ZK10
      DOUBLE PRECISION P(NPAR2DHP,NORD2DHP)
      DOUBLE PRECISION ANS1,ANS2,ANS3,ANS4,ANS5,ANS6,ANS7

      DOUBLE PRECISION SINZKZ(NORD2DHP+1),SSINZKZ(NORD2DHP+1)
      DOUBLE PRECISION COSZKZ(NORD2DHP+1),CCOSZKZ(NORD2DHP+1)

      COMPLEX*16 CZKZ(NORD2DHP+1)

      CHARACTER(64) COMMENT,FILEIN

      DATA ICAL/0/
      DATA PI/3.141592653589793D0/


C--- INITIALIZATION

      LUNIN=LUN2DHFIT
      FILEIN=FILE2DHFIT

      IF (ICAL.EQ.0) THEN

         OPEN(UNIT=LUNIN,FILE=FILEIN,STATUS='OLD')

            READ(LUN2DHFIT,'(A64)')COMMENT
            READ(LUN2DHFIT,*)PERLEN2DH,PHASE2DH,XYZ2DH
            READ(LUN2DHFIT,*)X2DHMIN,X2DHMAX
            READ(LUN2DHFIT,*)Y2DHMIN,Y2DHMAX
            READ(LUN2DHFIT,*)Z2DHMIN,Z2DHMAX
            READ(LUN2DHFIT,*)NORD2DH,NPAR2DH

          NPARTOT=(NORD2DH+1)/2*NPAR2DH*2
          DO IPAR=1,NPARTOT
         READ(LUN2DHFIT,*)PAR2DH(IPAR)
          ENDDO

            WRITE(LUNGFO,*)'     SUBROUTINE BPOLY2DH:'
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'     Fitparameters read from file:'
            WRITE(LUNGFO,*)'     ',FILEIN
            WRITE(LUNGFO,*)

            WRITE(LUNGFO,'(A64)')COMMENT
            WRITE(LUNGFO,*)
     &'     periodlength, phase, and scaling factor:'
            WRITE(LUNGFO,*)
     &'     (PERLEN2DH,PHASE2DH,XYZ2DH)'
            WRITE(LUNGFO,*)'     ',PERLEN2DH
            WRITE(LUNGFO,*)'     ',PHASE2DH
            WRITE(LUNGFO,*)'     ',XYZ2DH
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'     X2DHMIN,X2DHMAX:',X2DHMIN,X2DHMAX
            WRITE(LUNGFO,*)'     Y2DHMIN,Y2DHMAX:',Y2DHMIN,Y2DHMAX
            WRITE(LUNGFO,*)'     Z2DHMIN,Z2DHMAX:',Z2DHMIN,Z2DHMAX
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'     NORD2DH:',NORD2DH
            WRITE(LUNGFO,*)

          DO IPAR=1,NPARTOT
         WRITE(LUNGFO,*)PAR2DH(IPAR)
          ENDDO

          IF (NPAR2DH.GT.NPAR2DHP) THEN
         WRITE(LUNGFO,*)
     &'*** ERROR IN BPOLY2DH: DIMENSION NPAR2DHP IN BPOLY2DH.CMN EXCEEDED  ***'
         WRITE(LUNGFO,*)
         WRITE(6,*)
     &'*** ERROR IN BPOLY2DH: DIMENSION NPAR2DHP IN BPOLY2DH.CMN EXCEEDED  ***'
         WRITE(6,*)
         WRITE(6,*) '--- PROGRAM ABORTED DUE TO ERROR ---'
         STOP
          ENDIF

          IF (NORD2DH.GT.NORD2DHP) THEN
         WRITE(LUNGFO,*)
     &'*** ERROR IN BPOLY2DH: DIMENSION NORD2DHP IN BPOLY2DH.CMN EXCEEDED  ***'
         WRITE(LUNGFO,*)
         WRITE(6,*)
     &'*** ERROR IN BPOLY2DH: DIMENSION NORD2DHP IN BPOLY2DH.CMN EXCEEDED  ***'
         WRITE(6,*)
         WRITE(6,*) '--- PROGRAM ABORTED DUE TO ERROR ---'
         STOP
          ENDIF

         CLOSE(LUNIN)

      XOFF=PHASE2DH*PERLEN2DH
      ZL=PERLEN2DH*XYZ2DH

         ZK=2.D0*PI/ZL
         ZK2=ZK*ZK
         ZK4=ZK2*ZK2
         ZK6=ZK4*ZK2
         ZK8=ZK6*ZK2
         ZK10=ZK8*ZK2

         AXOUT=0.D0
         AYOUT=0.D0
         AZOUT=0.D0

         ICAL=1

      ENDIF !ICAL

C --- CHANGE COORDINATE SYSTEMS

      XX=-ZIN*XYZ2DH
      YY=YIN*XYZ2DH
      Z =(XIN+XOFF)*XYZ2DH

      XX2=XX*XX
      XX3=XX2*XX
      XX4=XX3*XX
      XX5=XX4*XX
      XX6=XX5*XX
      XX7=XX6*XX
      XX8=XX7*XX
      XX9=XX8*XX
      XX10=XX9*XX

      YY2=YY*YY
      YY3=YY2*YY
      YY4=YY3*YY
      YY5=YY4*YY
      YY6=YY5*YY
      YY7=YY6*YY
      YY8=YY7*YY
      YY9=YY8*YY
      YY10=YY9*YY

C--- MAGNETIC FIELD

      CZKZ(1)=CDEXP(DCMPLX(0.D0,ZK*Z))
      CZKZ(2)=CZKZ(1)*CZKZ(1)
      DO N=3,NORD2DH,2
            CZKZ(N)=CZKZ(N-2)*CZKZ(2)
      ENDDO !NORD2H

      DO N=1,NORD2DH,2
          SSINZKZ(N)=DIMAG(CZKZ(N))
          CCOSZKZ(N)=DREAL(CZKZ(N))
      ENDDO


      BXOUT=0.D0
      BYOUT=0.D0
      BZOUT=0.D0

      IPAR=0
      DO IHV=1,2

          IF (IHV.EQ.1) THEN

         X=XX
         X2=XX2
         X3=XX3
         X4=XX4
         X5=XX5
         X6=XX6
         X7=XX7
         X8=XX8
         X9=XX9
         X10=XX10

         Y=YY
         Y2=YY2
         Y3=YY3
         Y4=YY4
         Y5=YY5
         Y6=YY6
         Y7=YY7
         Y8=YY8
         Y9=YY9
         Y10=YY10

         DO N=1,NORD2DH,2
            sinzkz(N)=SSINZKZ(N)
            coszkz(N)=CCOSZKZ(N)
         ENDDO !N
          ELSE

         Y=XX
         Y2=XX2
         Y3=XX3
         Y4=XX4
         Y5=XX5
         Y6=XX6
         Y7=XX7
         Y8=XX8
         Y9=XX9
         Y10=XX10

         X=YY
         X2=YY2
         X3=YY3
         X4=YY4
         X5=YY5
         X6=YY6
         X7=YY7
         X7=YY7
         X8=YY8
         X9=YY9
         X10=YY10

C ACHTUNG VORZEICHENWECHSEL IN PARAMETERN WG ABLEITUNGEN VON SIN - COS !?
         DO N=1,NORD2DH,2
            sinzkz(N)=CCOSZKZ(N)
            coszkz(N)=SSINZKZ(N)
         ENDDO !N
          ENDIF

      DO N=1,NORD2DH,2

         N2=N*N
         N4=N2*N2
         N6=N4*N2
         N8=N6*N2
         N10=N8*N2

      DO NN=1,NPAR2DHP
          IPAR=IPAR+1
          P(NN,N)=PAR2DH(IPAR)
      ENDDO

c         INCLUDE 'red:poly2dh_bfeld.for'
      ans6=-0.00208333333333d0*p(2,n)*ZK6*n6*y8+0.116666666667d0*p
     . (2,n)*ZK4*n4*x6-0.116666666667d0*p(2,n)*ZK4*n4*y6+
     . 3.5*p(2,n)*ZK2*n2*x4-3.5*p(2,n)*ZK2*n2*y4+42.0*p
     . (2,n)*x2-42.0*p(2,n)*y2-0.0000115740740741d0*p(1,n)*ZK
     . 10*n10*x10-0.00104166666667d0*p(1,n)*ZK8*n8*x8-
     . 0.0583333333333d0*p(1,n)*ZK6*n6*x6-1.75*p(1,n)*ZK4*n
     . 4*x4-21.0*p(1,n)*ZK2*n2*x2-42.0*p(1,n)
      ans5=p(6,n)*x10-45.0*p(6,n)*x8*y2+210.0*p(6,n)*x6
     . *y4-210.0*p(6,n)*x4*y6+45.0*p(6,n)*x2*y8-p(6,n
     . )*y10-0.1*p(5,n)*ZK2*n2*x10+3.0*p(5,n)*ZK2*n2*x
     . 8*y2-7.0*p(5,n)*ZK2*n2*x6*y4+1.5*p(5,n)*ZK2*n
     . 2*x2*y8-0.0666666666667d0*p(5,n)*ZK2*n2*y10-3.0*p(
     . 5,n)*x8+84.0*p(5,n)*x6*y2-210.0*p(5,n)*x4*y4+
     . 84.0*p(5,n)*x2*y6-3.0*p(5,n)*y8+0.000555555555556d0*p(
     . 3,n)*ZK6*n6*x10-0.00625*p(3,n)*ZK6*n6*x8*y2+
     . 0.000138888888889d0*p(3,n)*ZK6*n6*y10+0.0375*p(3,n)*ZK
     . 4*n4*x8-0.35*p(3,n)*ZK4*n4*x6*y2+0.0125*p(3,n)*
     . ZK4*n4*y8+1.4*p(3,n)*ZK2*n2*x6-10.5*p(3,n)*ZK2*
     . n2*x4*y2+0.7*p(3,n)*ZK2*n2*y6+21.0*p(3,n)*x4-
     . 126.0*p(3,n)*x2*y2+21.0*p(3,n)*y4+0.0000231481481481d0
     . *p(2,n)*ZK8*n8*x10-0.0000231481481481d0*p(2,n)*ZK8*n
     . 8*y10+0.00208333333333d0*p(2,n)*ZK6*n6*x8+ans6
      ans4=0.0238095238095d0*ans5
      ans3=-ans4
      ans2=0.00000027557319224d0*(25920.0*(ZK2*n2*x2+ZK2*n2*y
     . 2+28.0)*(x2+4.0*x*y+y2)*(x2-4.0*x*y+y2)+432.0*(x8-
     . 14.0*x6*y2-14.0*x4*y4-14.0*x2*y6+y8)*ZK4*n4)*
     . (x+y)*(x-y)*p(4,n)+ans3
      ans7=COSZKz(n)
      ans1=ans2*ans7
      bxx=-ans1
      byy=(0.015873015873d0*(ZK2*n2*x2+ZK2*n2*y2+36.0)*(x2+
     . 2.0*x*y-y2)*(x2-2.0*x*y-y2)*(x+y)*(x-y)*p(5,n)-
     . 0.238095238095d0*(x4-2.0*x2*y2+0.2*y4)*(x4-10.0*x2*y
     . 2+5.0*y4)*p(6,n)-0.0000055114638448d0*(ZK8*n8*y8+72.0
     . *ZK6*n6*y6+3024.0*ZK4*n4*y4+60480.0*ZK2*n2*y2+
     . 3.6288e+5)*p(2,n)-0.0000055114638448d0*(432.0*(ZK4*n4*x4+
     . ZK4*n4*x2*y2+ZK4*n4*y4+42.0*ZK2*n2*x2+42.0*ZK
     . 2*n2*y2+840.0)+6.0*(x4+y4)*(x2+y2)*ZK6*n6)*(x
     . +y)*(x-y)*p(3,n)+0.0000055114638448d0*(5184.0*(x6-7.0*x2*y
     . 4+2.0*y6)*ZK2*n2+2.17728e+5*(x2-0.333333333333d0*y2)*
     . (x2-3.0*y2)+72.0*(x8-12.0*x2*y6+3.0*y8)*ZK4*n4)
     . *p(4,n))*COSZKz(n)*x*y
      ans5=231.0*p(3,n)*y4+0.0000231481481481d0*p(2,n)*ZK8*n8*
     . x10-0.00025462962963d0*p(2,n)*ZK8*n8*y10+
     . 0.0025462962963d0*p(2,n)*ZK6*n6*x8-0.0229166666667d0*p(2,
     . n)*ZK6*n6*y8+0.183333333333d0*p(2,n)*ZK4*n4*x6-
     . 1.28333333333d0*p(2,n)*ZK4*n4*y6+7.7*p(2,n)*ZK2*n2*x
     . 4-38.5*p(2,n)*ZK2*n2*y4+154.0*p(2,n)*x2-462.0*p(
     . 2,n)*y2-0.0000115740740741d0*p(1,n)*ZK10*n10*x10-
     . 0.00127314814815d0*p(1,n)*ZK8*n8*x8-0.0916666666667d0*
     &  p(1,n)*ZK6*n6*x6-3.85*p(1,n)*ZK4*n4*x4-77.0*p(1,n)*
     . ZK2*n2*x2-462.0*p(1,n)
      ans4=p(6,n)*x10-55.0*p(6,n)*x8*y2+330.0*p(6,n)*x6
     . *y4-462.0*p(6,n)*x4*y6+165.0*p(6,n)*x2*y8-11.0*p
     . (6,n)*y10-0.005*p(4,n)*ZK4*n4*x10+0.0916666666667d0*p
     . (4,n)*ZK4*n4*x8*y2-0.275*p(4,n)*ZK4*n4*x2*y8
     . +0.055*p(4,n)*ZK4*n4*y10-0.366666666667d0*p(4,n)*ZK2*n
     . 2*x8+6.6*p(4,n)*ZK2*n2*x6*y2-15.4*p(4,n)*ZK2*n
     . 2*x2*y6+3.3*p(4,n)*ZK2*n2*y8-13.2*p(4,n)*x6+
     . 277.2*p(4,n)*x4*y2-462.0*p(4,n)*x2*y4+92.4*p(4,n
     . )*y6+0.000555555555556d0*p(3,n)*ZK6*n6*x10-
     . 0.00763888888889d0*p(3,n)*ZK6*n6*x8*y2+0.00152777777778d0
     . *p(3,n)*ZK6*n6*y10+0.0458333333333d0*p(3,n)*ZK4*n4*x
     . 8-0.55*p(3,n)*ZK4*n4*x6*y2+0.1375*p(3,n)*ZK4*n
     . 4*y8+2.2*p(3,n)*ZK2*n2*x6-23.1*p(3,n)*ZK2*n2*x
     . 4*y2+7.7*p(3,n)*ZK2*n2*y6+46.2*p(3,n)*x4-462.0*p(
     . 3,n)*x2*y2+ans5
      ans3=0.0021645021645d0*ans4
      ans2=-ans3
      ans1=0.0000000250521083854d0*(8640.0*(x10-36.6666666667d0*x8*y
     . 2+110.0*x6*y4-55.0*x2*y8+7.33333333333d0*y10)*ZK2*n
     . 2+3.168e+5*(x6-33.0*x4*y2+27.0*x2*y4-3.0*y6)*(x
     . 2-3.0*y2))*p(5,n)+ans2
      ans6=SINZKz(n)*ZK*n*x
      bzz=ans1*ans6

C --- CHANGE COORDINATE SYSTEMS

          IF (IHV.EQ.1) THEN
             BXOUT=BXOUT+BZZ
             BYOUT=BYOUT+BYY
             BZOUT=BZOUT-BXX
          ELSE
             BZZ=-BZZ
             BXOUT=BXOUT+BZZ
             BYOUT=BYOUT+BXX
             BZOUT=BZOUT-BYY
          ENDIF

      ENDDO !NORD2DH
      ENDDO !IHV

      RETURN
      END
