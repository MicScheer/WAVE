*CMZ :  2.65/02 28/09/2009  08.45.56  by  Michael Scheer
*CMZ :  2.63/04 11/06/2009  12.17.42  by  Michael Scheer
*CMZ :  2.48/04 15/03/2004  16.23.51  by  Michael Scheer
*CMZ :  2.47/23 17/02/2004  12.09.11  by  Michael Scheer
*CMZ :  2.47/18 26/11/2003  18.23.46  by  Michael Scheer
*CMZ :  2.47/08 08/05/2003  12.49.51  by  Michael Scheer
*CMZ :  2.41/08 02/08/2002  19.47.58  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  14.40.26  by  Michael Scheer
*CMZ :  2.15/00 04/05/2000  17.01.34  by  Michael Scheer
*CMZ :  2.14/02 19/04/2000  17.02.46  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  11.45.40  by  Michael Scheer
*CMZ :  1.03/06 07/08/98  10.07.52  by  Michael Scheer
*CMZ : 00.02/04 12/02/97  11.01.21  by  Michael Scheer
*CMZ : 00.02/03 04/02/97  16.36.40  by  Michael Scheer
*-- Author :    Michael Scheer   23/01/97
        SUBROUTINE UTIL_LINEAR_FIT_USER(NARG,NFUN,NPAR,CURDAT,T)
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

      IMPLICIT NONE

*KEEP,bpoly2dh.
      include 'bpoly2dh.cmn'
*KEND.

      INTEGER ICAL
      INTEGER NARG,NFUN,NPAR,IPAR,IHV,NN
      INTEGER N
      DOUBLE PRECISION N2,N4,N6,N8,N10
      DOUBLE PRECISION CURDAT,T
      DIMENSION CURDAT(NARG+NFUN),T(NFUN,NPAR)

      DOUBLE PRECISION y,y2,y3,y4,y5,y6,Y7,Y8,Y9,Y10
      DOUBLE PRECISION YY,YY2,YY3,YY4,YY5,YY6,YY7,YY8,YY9,YY10
      DOUBLE PRECISION x,x2,x3,x4,x5,x6,X7,X8,X9,X10
      DOUBLE PRECISION XX,XX2,XX3,XX4,XX5,XX6,XX7,XX8,XX9,XX10
      DOUBLE PRECISION zl,zk,zkz,zk2,ZK4,ZK6,ZK8,ZK10
      DOUBLE PRECISION z,PI,ZOLD,XOLD,YOLD
      DOUBLE PRECISION COSZKZ(NORD2DHP+1), SINZKZ(NORD2DHP+1)
     &                  ,CCOSZKZ(NORD2DHP+1),SSINZKZ(NORD2DHP+1)

        COMPLEX*16 CZKZ(NORD2DHP+1)

      DATA ICAL/0/
      DATA PI/3.141592653589793D0/

      IF (ICAL.EQ.0) THEN
         ZL=PERLEN2DH
         zk=2.d0*pi/zl
         zk2=zk*zk
         ZK4=ZK2*ZK2
         ZK6=ZK4*ZK2
         ZK8=ZK6*ZK2
         ZK10=ZK8*ZK2

         XOLD=-1.D30
         YOLD=-1.D30
         ZOLD=-1.D30

         ICAL=1
      ENDIF !ICAL

      Xx=CURdat(1)
      Yy=CURdat(2)
       z=CURdat(3)

      IF (XX.NE.XOLD) THEN
         XX2=XX*XX
         XX3=XX2*XX
         XX4=XX3*XX
         XX5=XX4*XX
         XX6=XX5*XX
         XX7=XX6*XX
         XX8=XX7*XX
         XX9=XX8*XX
         XX10=XX9*XX
         XOLD=XX
      ENDIF

      IF (YY.NE.YOLD) THEN
         YY2=YY*YY
         YY3=YY2*YY
         YY4=YY3*YY
         YY5=YY4*YY
         YY6=YY5*YY
         YY7=YY6*YY
         YY8=YY7*YY
         YY9=YY8*YY
         YY10=YY9*YY
         YOLD=YY
      ENDIF

      IF (Z.NE.ZOLD) THEN
                zkz=zk*z
                CZKZ(1)=CDEXP(DCMPLX(0.D0,ZKZ))
                CZKZ(2)=CZKZ(1)*CZKZ(1)
         Ssinzkz(1)=DIMAG(CZKZ(1))
         Ccoszkz(1)=DREAL(CZKZ(1))
         DO N=3,NORD2DH,2
            CZKZ(N)=CZKZ(N-2)*CZKZ(2)
            Ssinzkz(N)=DIMAG(CZKZ(N))
            Ccoszkz(N)=DREAL(CZKZ(N))
         ENDDO !N
         ZOLD=Z
      ENDIF

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

c            INCLUDE 'red:poly2dh.for'
      tt(1,1)=-0.00000027557319224*(ZK10*n10*x10+90.0*ZK8*n
     . 8*x8+5040.0*ZK6*n6*x6+1.512e+5*ZK4*n4*x4+1.8144e+6
     . *ZK2*n2*x2+3.6288e+6)*COSZKz(n)
      tt(1,2)=0.00000055114638448*(ZK8*n8*x8+ZK8*n8*x6*y
     . 2+ZK8*n8*x4*y4+ZK8*n8*x2*y6+ZK8*n8*y8+
     . 90.0*ZK6*n6*x6+90.0*ZK6*n6*x4*y2+90.0*ZK6*n6*x
     . 2*y4+90.0*ZK6*n6*y6+5040.0*ZK4*n4*x4+5040.0*ZK4
     . *n4*x2*y2+5040.0*ZK4*n4*y4+1.512e+5*ZK2*n2*x2+
     . 1.512e+5*ZK2*n2*y2+1.8144e+6)*(x+y)*(x-y)*COSZKz(n)
      tt(1,3)=0.0000132275132275*(ZK6*n6*x10-11.25*ZK6*n6*x
     . 8*y2+0.25*ZK6*n6*y10+67.5*ZK4*n4*x8-630.0*ZK4*n
     . 4*x6*y2+22.5*ZK4*n4*y8+2520.0*ZK2*n2*x6-
     . 18900.0*ZK2*n2*x4*y2+1260.0*ZK2*n2*y6+37800.0*x4
     . -2.268e+5*x2*y2+37800.0*y4)*COSZKz(n)
      tt(1,4)=-0.00000027557319224*(25920.0*(ZK2*n2*x2+ZK2*n
     . 2*y2+28)*(x2+4.0*x*y+y2)*(x2-4.0*x*y+y2)+432.0*(
     . x8-14.0*x6*y2-14.0*x4*y4-14.0*x2*y6+y8)*ZK4*n
     . 4)*(x+y)*(x-y)*COSZKz(n)
      tt(1,5)=-0.00238095238095*(ZK2*n2*x10-30.0*ZK2*n2*x
     . 8*y2+70.0*ZK2*n2*x6*y4-15.0*ZK2*n2*x2*y8+
     . 0.666666666667*ZK2*n2*y10+30.0*x8-840.0*x6*y2+
     . 2100.0*x4*y4-840.0*x2*y6+30.0*y8)*COSZKz(n)
      tt(1,6)=0.0238095238095*(x4+4.0*x3*y-14.0*x2*y2+4.0*x
     . *y3+y4)*(x4-4.0*x3*y-14.0*x2*y2-4.0*x*y3+y4)*(
     . x+y)*(x-y)*COSZKz(n)
      tt(2,1)=0.0
      tt(2,2)=-0.0000055114638448*(ZK8*n8*y8+72.0*ZK6*n6*y
     . 6+3024.0*ZK4*n4*y4+60480.0*ZK2*n2*y2+3.6288e+5)*
     . COSZKz(n)*x*y
      tt(2,3)=-0.0000055114638448*(432.0*(ZK4*n4*x4+ZK4*n4*
     . x2*y2+ZK4*n4*y4+42.0*ZK2*n2*x2+42.0*ZK2*n2*y
     . 2+840)+6.0*(x4+y4)*(x2+y2)*ZK6*n6)*(x+y)*(x-y)*
     . COSZKz(n)*x*y
      tt(2,4)=0.0000055114638448*(5184.0*(x6-7.0*x2*y4+2.0*y
     . 6)*ZK2*n2+2.17728e+5*(x2-0.333333333333*y2)*(x2-3.0
     . *y2)+72.0*(x8-12.0*x2*y6+3.0*y8)*ZK4*n4)*COSZKz(n
     . )*x*y
      tt(2,5)=0.015873015873*(ZK2*n2*x2+ZK2*n2*y2+36)*(
     . x2+2.0*x*y-y2)*(x2-2.0*x*y-y2)*(x+y)*(x-y)*COSZKz(n)*x
     . *y
      tt(2,6)=-0.238095238095*(x4-2.0*x2*y2+0.2*y4)*(x4-
     . 10.0*x2*y2+5.0*y4)*COSZKz(n)*x*y
      tt(3,1)=0.0000000250521083854*(ZK10*n10*x10+110.0*ZK8*n
     . 8*x8+7920.0*ZK6*n6*x6+3.3264e+5*ZK4*n4*x4+
     . 6.6528e+6*ZK2*n2*x2+3.99168e+7)*SINZKz(n)*ZK*n*x
      tt(3,2)=-0.0000000501042167709*(ZK8*n8*x10-11.0*ZK8*n
     . 8*y10+110.0*ZK6*n6*x8-990.0*ZK6*n6*y8+7920.0*ZK4*
     . n4*x6-55440.0*ZK4*n4*y6+3.3264e+5*ZK2*n2*x4-
     . 1.6632e+6*ZK2*n2*y4+6.6528e+6*x2-1.99584e+7*y2)*
     . SINZKz(n)*ZK*n*x
      tt(3,3)=-0.0000012025012025*(ZK6*n6*x10-13.75*ZK6*n6*
     . x8*y2+2.75*ZK6*n6*y10+82.5*ZK4*n4*x8-990.0*ZK4*
     . n4*x6*y2+247.5*ZK4*n4*y8+3960.0*ZK2*n2*x6-
     . 41580.0*ZK2*n2*x4*y2+13860.0*ZK2*n2*y6+83160.0*x
     . 4-8.316e+5*x2*y2+4.158e+5*y4)*SINZKz(n)*ZK*n*x
      tt(3,4)=0.0000108225108225*(ZK4*n4*x10-18.3333333333*ZK
     . 4*n4*x8*y2+55.0*ZK4*n4*x2*y8-11.0*ZK4*n4*y10
     . +73.3333333333*ZK2*n2*x8-1320.0*ZK2*n2*x6*y2+
     . 3080.0*ZK2*n2*x2*y6-660.0*ZK2*n2*y8+2640.0*x6-
     . 55440.0*x4*y2+92400.0*x2*y4-18480.0*y6)*SINZKz(n)*ZK*
     . n*x
      tt(3,5)=0.0000000250521083854*(8640.0*(x10-36.6666666667*x
     . 8*y2+110.0*x6*y4-55.0*x2*y8+7.33333333333*y10)*ZK
     . 2*n2+3.168e+5*(x6-33.0*x4*y2+27.0*x2*y4-3.0*y6
     . )*(x2-3.0*y2))*SINZKz(n)*ZK*n*x
      tt(3,6)=-0.0021645021645*(x10-55.0*x8*y2+330.0*x6*y
     . 4-462.0*x4*y6+165.0*x2*y8-11.0*y10)*SINZKz(n)*ZK*n*x


          IPAR=IPAR+1
          IF (IHV.EQ.1) THEN
             T(1,IPAR)=TT(1,NN)
             T(2,IPAR)=TT(2,NN)
             T(3,IPAR)=TT(3,NN)
          ELSE
             T(1,IPAR)=TT(2,NN)
             T(2,IPAR)=TT(1,NN)
             T(3,IPAR)=-TT(3,NN)
          ENDIF

      ENDDO !NN
      ENDDO !N
      ENDDO !IHV


      RETURN
      END
