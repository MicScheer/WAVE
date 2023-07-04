*CMZ :  2.41/10 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  17.27.12  by  Michael Scheer
*CMZ :  2.14/02 19/04/2000  17.02.45  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  11.45.40  by  Michael Scheer
*CMZ : 00.02/01 17/12/96  12.01.36  by  Michael Scheer
*-- Author :    Michael Scheer   16/12/96

      SUBROUTINE BHELITRIESTE(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)

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
*KEND.

      INTEGER ICAL,LUNIN

      DOUBLE PRECISION PI
      DOUBLE PRECISION XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT
      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,BXH,BXV,BYH,BYV,BZH,BZV
      DOUBLE PRECISION X2,X4,X6,Y2,Y4,Y6
      DOUBLE PRECISION ZL,ZK,ZK2,SINZKZ,SINZKZ3,COSZKZ,COSZKZ3

      COMPLEX*16 CZKZ,CZKZ3

      DOUBLE PRECISION BH1,B2H1,C3H1,D4H1
      DOUBLE PRECISION BH3,B2H3,C3H3,D4H3
      DOUBLE PRECISION BV1,B2V1,C3V1,D4V1
      DOUBLE PRECISION BV3,B2V3,C3V3,D4V3

      DOUBLE PRECISION A2H1,B3H1,C4H1,A3H1,B4H1,A4H1
      DOUBLE PRECISION A2V1,B3V1,C4V1,A3V1,B4V1,A4V1
      DOUBLE PRECISION A2H3,B3H3,C4H3,A3H3,B4H3,A4H3
      DOUBLE PRECISION A2V3,B3V3,C4V3,A3V3,B4V3,A4V3

      CHARACTER(64) COMMENT,FILEIN

      DATA ICAL/0/

      DATA LUNIN/10/
      DATA FILEIN/'HELI_TRIESTE.FIT'/
      DATA PI/3.141592653589793D0/

C--- INITIALIZATION

      IF (ICAL.EQ.0) THEN

         OPEN(UNIT=LUNIN,FILE=FILEIN,STATUS='OLD')

            READ(LUNIN,'(A64)')COMMENT

            READ(LUNIN,*)ZL

            READ(LUNIN,*)BH1
            READ(LUNIN,*)B2H1
            READ(LUNIN,*)C3H1
            READ(LUNIN,*)D4H1
            READ(LUNIN,*)BH3
            READ(LUNIN,*)B2H3
            READ(LUNIN,*)C3H3
            READ(LUNIN,*)D4H3

            READ(LUNIN,*)BV1
            READ(LUNIN,*)B2V1
            READ(LUNIN,*)C3V1
            READ(LUNIN,*)D4V1
            READ(LUNIN,*)BV3
            READ(LUNIN,*)B2V3
            READ(LUNIN,*)C3V3
            READ(LUNIN,*)D4V3

            WRITE(LUNGFO,*)'     SUBROUTINE BHELITRIESTE:'
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'     COEFFICENT FILE:'
            WRITE(LUNGFO,'(''      '',A64)')FILEIN
            WRITE(LUNGFO,*)'     COMMENT:'
            WRITE(LUNGFO,'(''      '',A64)')COMMENT
            WRITE(LUNGFO,*)

            WRITE(LUNGFO,*)'     Z-Lambda:',ZL
            WRITE(LUNGFO,*)

            WRITE(LUNGFO,*)'     COEFFICIENTS:'
            WRITE(LUNGFO,*)'     ',BH1
            WRITE(LUNGFO,*)'     ',B2H1
            WRITE(LUNGFO,*)'     ',C3H1
            WRITE(LUNGFO,*)'     ',D4H1
            WRITE(LUNGFO,*)'     ',BH3
            WRITE(LUNGFO,*)'     ',B2H3
            WRITE(LUNGFO,*)'     ',C3H3
            WRITE(LUNGFO,*)'     ',D4H3

            WRITE(LUNGFO,*)'     ',BV1
            WRITE(LUNGFO,*)'     ',B2V1
            WRITE(LUNGFO,*)'     ',C3V1
            WRITE(LUNGFO,*)'     ',D4V1
            WRITE(LUNGFO,*)'     ',BV3
            WRITE(LUNGFO,*)'     ',B2V3
            WRITE(LUNGFO,*)'     ',C3V3
            WRITE(LUNGFO,*)'     ',D4V3

         CLOSE(LUNIN)

         ZK=2.D0*PI/ZL
         ZK2=ZK*ZK

         A2H1=(ZK2     - 2.D0*B2H1)/6.D0
         B3H1=(ZK2*B2H1-12.D0*C3H1)/6.D0
         C4H1=(ZK2*C3H1-30.D0*D4H1)/6.D0
         A3H1=(ZK2*A2H1- 2.D0*B3H1)/20.D0
         B4H1=(ZK2*B3H1-12.D0*C4H1)/20.D0
         A4H1=(ZK2*A3H1- 2.D0*B4H1)/42.D0

         A2H3=(9.D0*ZK2     - 2.D0*B2H3)/6.D0
         B3H3=(9.D0*ZK2*B2H3-12.D0*C3H3)/6.D0
         C4H3=(9.D0*ZK2*C3H3-30.D0*D4H3)/6.D0
         A3H3=(9.D0*ZK2*A2H3- 2.D0*B3H3)/20.D0
         B4H3=(9.D0*ZK2*B3H3-12.D0*C4H3)/20.D0
         A4H3=(9.D0*ZK2*A3H3- 2.D0*B4H3)/42.D0

         A2V1=(ZK2     - 2.D0*B2V1)/6.D0
         B3V1=(ZK2*B2V1-12.D0*C3V1)/6.D0
         C4V1=(ZK2*C3V1-30.D0*D4V1)/6.D0
         A3V1=(ZK2*A2V1- 2.D0*B3V1)/20.D0
         B4V1=(ZK2*B3V1-12.D0*C4V1)/20.D0
         A4V1=(ZK2*A3V1- 2.D0*B4V1)/42.D0

         A2V3=(9.D0*ZK2     - 2.D0*B2V3)/6.D0
         B3V3=(9.D0*ZK2*B2V3-12.D0*C3V3)/6.D0
         C4V3=(9.D0*ZK2*C3V3-30.D0*D4V3)/6.D0
         A3V3=(9.D0*ZK2*A2V3- 2.D0*B3V3)/20.D0
         B4V3=(9.D0*ZK2*B3V3-12.D0*C4V3)/20.D0
         A4V3=(9.D0*ZK2*A3V3- 2.D0*B4V3)/42.D0

         AXOUT=0.D0
         AYOUT=0.D0
         AZOUT=0.D0

         ICAL=1

      ENDIF !ICAL

C --- CHANGE COORDINATE SYSTEMS

      X=-ZIN
      Y=YIN
      Z=XIN

      X2=X*X
      X4=X2*X2
      X6=X4*X2

      Y2=Y*Y
      Y4=Y2*Y2
      Y6=Y4*Y2

C--- MAGNETIC FIELD

      CZKZ=CDEXP(DCMPLX(0.D0,ZK*Z))
      CZKZ3=CZKZ*CZKZ*CZKZ
      SINZKZ =DIMAG( CZKZ)
      SINZKZ3=DIMAG(CZKZ3)
      COSZKZ =DREAL( CZKZ)
      COSZKZ3=DREAL(CZKZ3)
C--------------------------------------------- s HELI_TRIESTE.RED,~.COM
      bxv=-2.D0*((2.D0*(c3v1+c4v1*y2)*x2+b2v1+(b3v1+b4v1*y2)*y2
     . )*SINZKZ*bv1+(2.D0*(c3v3+c4v3*y2)*x2+b2v3+(b3v3+b4v3*y
     . 2)*y2)*SINZKZ3*bv3)*x*y
      byv=-(3.D0*SINZKZ*a2v1*bv1*y2+5.D0*SINZKZ*a3v1*bv1*y4+
     . 7.D0*SINZKZ*a4v1*bv1*y6+SINZKZ*b2v1*bv1*x2+3.D0*SINZ
     . KZ*b3v1*bv1*x2*y2+5.D0*SINZKZ*b4v1*bv1*x2*y4+SINZ
     . KZ*bv1*c3v1*x4+3.D0*SINZKZ*bv1*c4v1*x4*y2+SINZKZ*
     . bv1*d4v1*x6+SINZKZ*bv1+3.D0*SINZKZ3*a2v3*bv3*y2+5.D0*
     . SINZKZ3*a3v3*bv3*y4+7.D0*SINZKZ3*a4v3*bv3*y6+
     . SINZKZ3*b2v3*bv3*x2+3.D0*SINZKZ3*b3v3*bv3*x2*y
     . 2+5.D0*SINZKZ3*b4v3*bv3*x2*y4+SINZKZ3*bv3*
     . c3v3*x4+3.D0*SINZKZ3*bv3*c4v3*x4*y2+SINZKZ3*
     . bv3*d4v3*x6+SINZKZ3*bv3)
      bzv=-((d4v1*x6+1.D0+c4v1*x4*y2+c3v1*x4+b4v1*x2*y4+b3v1
     . *x2*y2+b2v1*x2+a4v1*y6+a3v1*y4+a2v1*y2)*COSZKZ*
     . bv1+3.D0*(d4v3*x6+1.D0+c4v3*x4*y2+c3v3*x4+b4v3*x2*y4+
     . b3v3*x2*y2+b2v3*x2+a4v3*y6+a3v3*y4+a2v3*y2)*COS
     . ZKZ3*bv3)*y*zk
      bxh=-(3.D0*COSZKZ*a2h1*bh1*x2+5.D0*COSZKZ*a3h1*bh1*x4+
     . 7.D0*COSZKZ*a4h1*bh1*x6+COSZKZ*b2h1*bh1*y2+3.D0*SINZ
     . KZ*b3h1*bh1*x2*y2+5.D0*COSZKZ*b4h1*bh1*x4*y2+SINZ
     . KZ*bh1*c3h1*y4+3.D0*COSZKZ*bh1*c4h1*x2*y4+COSZKZ*
     . bh1*d4h1*y6+COSZKZ*bh1+3.D0*COSZKZ3*a2h3*bh3*x2+5.D0*
     . COSZKZ3*a3h3*bh3*x4+7.D0*COSZKZ3*a4h3*bh3*x6+
     . COSZKZ3*b2h3*bh3*y2+3.D0*COSZKZ3*b3h3*bh3*x2*y
     . 2+5.D0*COSZKZ3*b4h3*bh3*x4*y2+COSZKZ3*bh3*
     . c3h3*y4+3.D0*COSZKZ3*bh3*c4h3*x2*y4+COSZKZ3*
     . bh3*d4h3*y6+COSZKZ3*bh3)
      byh=-2.D0*((2.D0*(c3h1+c4h1*x2)*y2+b2h1+(b3h1+b4h1*x2)*x2
     . )*COSZKZ*bh1+(2.D0*(c3h3+c4h3*x2)*y2+b2h3+(b3h3+b4h3*x
     . 2)*x2)*COSZKZ3*bh3)*x*y
      bzh=((d4h1*y6+1.D0+c4h1*x2*y4+c3h1*y4+b4h1*x4*y2+b3h1*
     . x2*y2+b2h1*y2+a4h1*x6+a3h1*x4+a2h1*x2)*SINZKZ*
     . bh1+3.D0*(d4h3*y6+1.D0+c4h3*x2*y4+c3h3*y4+b4h3*x4*y2+
     . b3h3*x2*y2+b2h3*y2+a4h3*x6+a3h3*x4+a2h3*x2)*SIN
     . ZKZ3*bh3)*x*zk
C---------------------------------------------

      BX=BXH+BXV
      BY=BYH+BYV
      BZ=BZH+BZV

C --- CHANGE COORDINATE SYSTEMS

      BXOUT=BZ
      BYOUT=BY
      BZOUT=-BX

      RETURN
      END
