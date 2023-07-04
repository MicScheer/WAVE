*CMZ :  2.41/10 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.16/04 19/06/2000  14.36.05  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  17.24.47  by  Michael Scheer
*CMZ :  2.14/02 19/04/2000  17.02.45  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  11.45.40  by  Michael Scheer
*CMZ : 00.02/03 21/01/97  14.52.55  by  Michael Scheer
*CMZ : 00.02/01 17/12/96  12.01.36  by  Michael Scheer
*-- Author :    Michael Scheer   16/12/96

      SUBROUTINE BTOSCA(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)
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
C *** ATTENTION: HALBACH ANSATZ CAN NOT BE FITTED SINCE FIELD MUST BE
C         BY ~ SIN(K*Z)

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
      DOUBLE PRECISION X,Y,Z,BXX,BYY,BZZ
      DOUBLE PRECISION X2,X3,X4,X5,X6,Y2,Y3,Y4,Y5,Y6
      DOUBLE PRECISION ZL,ZK,ZK2,ZK22,ZK23,SIN3ZKZ,SINZKZ,COSZKZ,COS3ZKZ

      DOUBLE PRECISION ANS1,ANS2
     &,bv01
     &,b2v1
     &,c3v1
     &,d4v1
     &,bv03
     &,b2v3
     &,c3v3
     &,d4v3
     &,bH01
     &,b2H1
     &,c3H1
     &,d4H1
     &,bH03
     &,b2H3
     &,c3H3
     &,d4H3

      COMPLEX*16 CZKZ,CZKZ3

      CHARACTER(64) COMMENT,FILEIN

      DATA ICAL/0/

      DATA LUNIN/10/
      DATA FILEIN/'TOSCA.FIT'/
      DATA PI/3.141592653589793D0/

C--- INITIALIZATION

      IF (ICAL.EQ.0) THEN

         OPEN(UNIT=LUNIN,FILE=FILEIN,STATUS='OLD')

            READ(LUNIN,'(A64)')COMMENT

          READ(LUNIN,*)ZL

               READ(LUNIN,*)bH01
               READ(LUNIN,*)b2H1
               READ(LUNIN,*)c3H1
               READ(LUNIN,*)d4H1
               READ(LUNIN,*)bH03
               READ(LUNIN,*)b2H3
               READ(LUNIN,*)c3H3
               READ(LUNIN,*)d4H3

               READ(LUNIN,*)bv01
               READ(LUNIN,*)b2v1
               READ(LUNIN,*)c3v1
               READ(LUNIN,*)d4v1
               READ(LUNIN,*)bv03
               READ(LUNIN,*)b2v3
               READ(LUNIN,*)c3v3
               READ(LUNIN,*)d4v3

            WRITE(LUNGFO,*)'     SUBROUTINE BTOSCA:'
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'     COEFFICENT FILE:'
            WRITE(LUNGFO,'(''      '',A64)')FILEIN
            WRITE(LUNGFO,*)'     COMMENT:'
            WRITE(LUNGFO,'(''      '',A64)')COMMENT
            WRITE(LUNGFO,*)

            WRITE(LUNGFO,*)'     Z-Lambda:',ZL
            WRITE(LUNGFO,*)

            WRITE(LUNGFO,*)'     COEFFICIENTS:'

              WRITE(LUNGFO,*)bH01
              WRITE(LUNGFO,*)b2H1
              WRITE(LUNGFO,*)c3H1
              WRITE(LUNGFO,*)d4H1
              WRITE(LUNGFO,*)bH03
              WRITE(LUNGFO,*)b2H3
              WRITE(LUNGFO,*)c3H3
              WRITE(LUNGFO,*)d4H3

              WRITE(LUNGFO,*)bv01
              WRITE(LUNGFO,*)b2v1
              WRITE(LUNGFO,*)c3v1
              WRITE(LUNGFO,*)d4v1
              WRITE(LUNGFO,*)bv03
              WRITE(LUNGFO,*)b2v3
              WRITE(LUNGFO,*)c3v3
              WRITE(LUNGFO,*)d4v3

         CLOSE(LUNIN)

         ZK=2.D0*PI/ZL
         ZK2=ZK*ZK
         ZK22=ZK2*ZK2
         ZK23=ZK2*ZK2*ZK2

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
      X3=X2*X
      X4=X3*X
      X5=X4*X
      X6=X5*X

      Y2=Y*Y
      Y3=Y2*Y
      Y4=Y3*Y
      Y5=Y4*Y
      Y6=Y5*Y

C--- MAGNETIC FIELD


      CZKZ=CDEXP(DCMPLX(0.D0,ZK*Z))
      CZKZ3=CZKZ*CZKZ*CZKZ

      SINZKZ =DIMAG( CZKZ)
      SIN3ZKZ=DIMAG(CZKZ3)

      COSZKZ =DREAL( CZKZ)
      COS3ZKZ=DREAL(CZKZ3)


C------------------ s TRIESTE_TOSCA_... IN RED:REDUCE.CMZ
C     include 'red:trieste_tosca_bfeld.for'
C     ERZEUGEN MIT TRIESTE_TOSCA_RED AUS REDUCE.CMZ
C       BEARBEITEN MIT TRIESTE_TOSCA_EDI AUS REDUCE.CMZ
C
      ans2=0.00138888888889*(648.0*(x6*zk2+1.11111111111*y4)+
     . 3240.0*(y2*zk2-1.33333333333)*x2*y2-6480.0*(y2*zk2-
     . 0.111111111111)*x4)*c3h3*cos3zkz-0.00138888888889*(6.0*(x
     . 6*zk22-120.0*y2)-30.0*(y2*zk2-4.0)*x4*zk2-360.0*(y2*
     . zk2-2.0)*x2)*b2h1*coszkz-0.00138888888889*(486.0*(x6*zk2
     . 2-1.48148148148*y2)-2430.0*(y2*zk2-0.444444444444)*x4*
     . zk2-3240.0*(y2*zk2-0.222222222222)*x2)*b2h3*cos3zkz
      ans1=(0.0166666666667*(b2v1*sinzkz*y4*zk22+20.0*b2v1*
     . sinzkz*y2*zk2+120.0*b2v1*sinzkz+81.0*b2v3*sin3zkz*y4*
     . zk22+180.0*b2v3*sin3zkz*y2*zk2+120.0*b2v3*sin3zkz)+6.0*(
     . d4v1*sinzkz+d4v3*sin3zkz)*(x2-0.333333333333*y2)*(x2-
     . 3.0*y2))*x*y-(d4h1*coszkz+d4h3*cos3zkz)*(x2+4.0*x*y+y2
     . )*(x2-4.0*x*y+y2)*(x+y)*(x-y)+0.00138888888889*(x6*zk2
     . 3+30.0*x4*zk22+360.0*x2*zk2+720.0)*bh01*coszkz+1.0125*(
     . x6*zk23+3.33333333333*x4*zk22+4.44444444444*x2*zk2+
     . 0.987654320988)*bh03*cos3zkz-0.00138888888889*(288.0*(y2*
     . zk2+10.0)*y2-480.0*(y2*zk2+6.0)*x2)*c3v1*sinzkz*x*y-
     . 0.00138888888889*(2592.0*(y2*zk2+1.11111111111)*y2-4320.0*
     . (y2*zk2+0.666666666667)*x2)*c3v3*sin3zkz*x*y+
     . 0.00138888888889*(72.0*(x6*zk2+10.0*y4)+360.0*(y2*zk2-
     . 12.0)*x2*y2-720.0*(y2*zk2-1.0)*x4)*c3h1*coszkz+ans2
      bxx=-ans1
      ans2=-0.00138888888889*(486.0*(y4*zk22+2.22222222222*y2*
     . zk2+1.48148148148)*y2-2430.0*(y4*zk22+1.33333333333*y2
     . *zk2+0.296296296296)*x2)*b2v3*sin3zkz+0.00138888888889*(
     . 72.0*(y2*zk2+10.0)*y4+360.0*(y2*zk2+2.0)*x4-720.0*(y
     . 2*zk2+6.0)*x2*y2)*c3v1*sinzkz+0.00138888888889*(648.0*(y
     . 2*zk2+1.11111111111)*y4+3240.0*(y2*zk2+0.222222222222)*x
     . 4-6480.0*(y2*zk2+0.666666666667)*x2*y2)*c3v3*sin3zkz
      ans1=0.00138888888889*(bv01*sinzkz*y6*zk23+30.0*bv01*
     . sinzkz*y4*zk22+360.0*bv01*sinzkz*y2*zk2+720.0*bv01*
     . sinzkz+729.0*bv03*sin3zkz*y6*zk23+2430.0*bv03*sin3zkz*y
     . 4*zk22+3240.0*bv03*sin3zkz*y2*zk2+720.0*bv03*sin3zkz)+
     . 6.0*(d4h1*coszkz+d4h3*cos3zkz)*(x2-0.333333333333*y2)*(x
     . 2-3.0*y2)*x*y+(d4v1*sinzkz+d4v3*sin3zkz)*(x2+4.0*x*y+y
     . 2)*(x2-4.0*x*y+y2)*(x+y)*(x-y)-0.4*(x4*zk2-
     . 1.66666666667*x2*y2*zk2+10.0*x2-10.0*y2)*c3h1*coszkz*
     . x*y-3.6*(x4*zk2-1.66666666667*x2*y2*zk2+1.11111111111*x
     . 2-1.11111111111*y2)*c3h3*cos3zkz*x*y+0.0166666666667*(x
     . 4*zk22+20.0*x2*zk2+120.0)*b2h1*coszkz*x*y+1.35*(x4*zk2
     . 2+2.22222222222*x2*zk2+1.48148148148)*b2h3*cos3zkz*x*y-
     . 0.00138888888889*(6.0*(y4*zk22+20.0*y2*zk2+120.0)*y2-
     . 30.0*(y4*zk22+12.0*y2*zk2+24.0)*x2)*b2v1*sinzkz+ans2
      byy=-ans1
      ans2=((0.2*(y2*zk2+10.0)*y2-0.166666666667*(y2*zk2+6.0)*x
     . 2)*x2-0.0142857142857*(y2*zk2+14.0)*y4)*c3v1*coszkz*y
     . +((5.4*(y2*zk2+1.11111111111)*y2-4.5*(y2*zk2+
     . 0.666666666667)*x2)*x2-0.385714285714*(y2*zk2+
     . 1.55555555556)*y4)*c3v3*cos3zkz*y+(0.385714285714*(x6*zk2
     . +7.77777777778*y4)+4.5*(y2*zk2-1.33333333333)*x2*y2-
     . 5.4*(y2*zk2-0.111111111111)*x4)*c3h3*sin3zkz*x+(
     . 0.0142857142857*x6*zk2+y4+0.166666666667*(y2*zk2-12.0)*x
     . 2*y2-0.2*(y2*zk2-1.0)*x4)*c3h1*sinzkz*x-(
     . 0.00119047619048*(x6*zk22-840.0*y2)-0.00833333333333*(y
     . 2*zk2-4.0)*x4*zk2-0.166666666667*(y2*zk2-2.0)*x2)*b2h1
     . *sinzkz*x-(0.289285714286*(x6*zk22-10.3703703704*y2)-
     . 2.025*(y2*zk2-0.444444444444)*x4*zk2-4.5*(y2*zk2-
     . 0.222222222222)*x2)*b2h3*sin3zkz*x
      ans1=0.000198412698413*(bh01*sinzkz*x6*zk23+42.0*bh01*
     . sinzkz*x4*zk22+840.0*bh01*sinzkz*x2*zk2+5040.0*bh01*
     . sinzkz+2187.0*bh03*sin3zkz*x6*zk23+10206.0*bh03*sin3zkz
     . *x4*zk22+22680.0*bh03*sin3zkz*x2*zk2+15120.0*bh03*
     . sin3zkz)*x-(d4v1*coszkz+3.0*d4v3*cos3zkz)*(x6-5.0*x4*y
     . 2+3.0*x2*y4-0.142857142857*y6)*y-0.142857142857*(d4h1*
     . sinzkz+3.0*d4h3*sin3zkz)*(x6-21.0*x4*y2+35.0*x2*y4
     . -7.0*y6)*x-0.000198412698413*(bv01*coszkz*y6*zk23+42.0*
     . bv01*coszkz*y4*zk22+840.0*bv01*coszkz*y2*zk2+5040.0*
     . bv01*coszkz+2187.0*bv03*cos3zkz*y6*zk23+10206.0*bv03*
     . cos3zkz*y4*zk22+22680.0*bv03*cos3zkz*y2*zk2+15120.0*
     . bv03*cos3zkz)*y+(0.00119047619048*(y4*zk22+28.0*y2*zk2+
     . 280.0)*y2-0.00833333333333*(y4*zk22+20.0*y2*zk2+120.0)
     . *x2)*b2v1*coszkz*y+(0.289285714286*(y4*zk22+
     . 3.11111111111*y2*zk2+3.45679012346)*y2-2.025*(y4*zk22+
     . 2.22222222222*y2*zk2+1.48148148148)*x2)*b2v3*cos3zkz*y+
     . ans2
      bzz=ans1*zk
C
C---------------------------------------------

C --- CHANGE COORDINATE SYSTEMS

      BXOUT=BZZ
      BYOUT=BYY
      BZOUT=-BXX

      RETURN
      END
