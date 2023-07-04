*CMZ :  3.00/00 11/03/2013  10.38.00  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.61/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.41.44  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.01/10 03/09/96  15.06.58  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  16.24.37  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.50.18  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.18  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE ERZFUN  (GAMMA,ZI,BXI,BYI,BZI,BXF,BYF,BZF,
     &                     AXI,AYI,AZI,AXF,AYF,AZF,XI,XPI,YI,YPI,
     &            XF,XPF,YF,YPF,X0,XF0)
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


c     The subroutine calculates the phase space coordinates according
c     to the coefficients of the generating function on the data file
c     WAVE_ERZFUN.IN. The indices I,J,K,L of the coefficients A(I,J,K,L)
c     are considered as mathematical indices rather than FORTRAN
c     array indices.

c     INPUT: XI,XPI,YI,YPI,AXI,AYI
c     OUTPUT:XI,XPI,YI,YPI,XF,XPF,YF,YPF,AXI,AYI;
c              i.e. the input is overwritten !!

c     INPUT-FILE: WAVE_ERZFUN.IN
c     LUN=LIN


      IMPLICIT NONE

      CHARACTER(65) TRANCODE

      INTEGER LIN,ICAL,I,J,K,L,IDUM,nkoef

*KEEP,genfun.
      include 'genfun.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DOUBLE PRECISION A(NORDNG+1,NORDNG+1,NORDNG+1,NORDNG+1),PXF,PYF,PXI,PYI,
     &         BRHO,PEL,GAMMA
     &        ,ZI,AXI,AYI,AXF,AYF,BXI,BYI,BXF,BYF
     &        , XI,XPI,YI,YPI,XF,XPF,YF,YPF,XII,PXFJ,YIK,PYFL,DSQ1,
     &          BZI,AZI,AZF,BZF,AX0,AY0,AZ0,AXF0,AYF0,AZF0,X0,XF0
     &        ,BX0,BY0,BZ0,BXF0,BYF0,BZF0,DUM

      DATA LIN/20/
      DATA ICAL/0/

      IF (ICAL.NE.-1) THEN

          PEL=EMASSE1*DSQRT( (GAMMA+1.D0)*(GAMMA-1.D0) )
          BRHO=PEL/CLIGHT1 !ABSOLUTE VALUE

          DO I=0,NORDNG
          DO J=0,NORDNG
          DO K=0,NORDNG
          DO L=0,NORDNG

         A(I+1,J+1,K+1,L+1)=0.D0

          ENDDO
          ENDDO
          ENDDO
          ENDDO

          OPEN(UNIT=LIN,FILE='wave_erzfun.in',FORM='FORMATTED',STATUS='OLD')

          READ(LIN,'(1A65)')TRANCODE
          read(lin,*)idum
          read(lin,*)dum
          read(lin,*)dum
          read(lin,*)dum
          read(lin,*)dum
          read(lin,*)dum
          read(lin,*)dum
          read(lin,*)dum
          read(lin,*)dum
          read(lin,*)dum
          read(lin,*)dum
          read(lin,*)dum
          read(lin,*)dum

          DO WHILE (.TRUE.)
               READ (LIN,*,END=99) I,J,K,L,A(I+1,J+1,K+1,L+1)
          END DO

99        CLOSE(LIN)

          IF ( A(2,1,1,1).NE.0. .OR. A(1,2,1,1).NE.0. .OR. A(1,1,2,1)
     &           .NE.0. .OR. A(1,1,1,2).NE.0.) THEN

           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)
     &       '     ***  ERZFUN: CLOSED ORBIT SET TO ZERO ***'
           WRITE(LUNGFO,*)
           WRITE(6,*)
           WRITE(6,*) '***  ERZFUN: CLOSED ORBIT SET TO ZERO ***'
           WRITE(6,*)

           A(2,1,1,1)=0.
           A(1,2,1,1)=0.
           A(1,1,2,1)=0.
           A(1,1,1,2)=0.

           ENDIF

           CALL MYBFELD(X0,0.D0,0.D0,BXF0,BYF0,BZF0,AXF0,AYF0,AZF0)
           CALL MYBFELD(XF0,0.D0,0.D0,BX0,BY0,BZ0,AX0,AY0,AZ0)

           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'      ERZFUN called: Code on file wave_erzfun.in:'
           WRITE(LUNGFO,*)'      ',TRANCODE
           WRITE(LUNGFO,*)

           ICAL=-1

         ENDIF

C--- NON-CANONICAL VARIABLES, THE

      XPF=XPI
      YPF=YPI

C--- CANONICAL VARIABLES

      BXF=BXI
      BYF=BYI
      BZF=BZI

      AXF=AXI-AX0
      AYF=AYI-AY0
      AZF=AZI-AZ0

      DSQ1=1.D0/DSQRT(1.D0+XPI*XPI+YPI*YPI)

C--- THE VECTOR POTENTIAL IS MORE OR LESS ARBITRARY BUT PXF,AZF AND XPF
C    MUST BE CONSISTENT

      PXF=-AZF/BRHO+XPI*DSQ1
      PYF=-AYF/BRHO+YPI*DSQ1

      XF=0.D0
      DO I=0,NORDNG
          DO J=1,NORDNG
         DO K=0,NORDNG
             DO L=0,NORDNG

      IF( I+J+K+L .LE. NORDNG ) THEN

         IF(I.EQ.0 .AND. XI.EQ. 0.D0)  THEN
            XII=1.D0
         ELSE
            XII=XI**I
         ENDIF
         IF(J-1.EQ.0 .AND. PXF.EQ. 0.D0)  THEN
            PXFJ=1.D0
         ELSE
            PXFJ=PXF**(J-1)
         ENDIF
         IF(K.EQ.0 .AND. YI.EQ.0.D0) THEN
            YIK=1.D0
         ELSE
            YIK=YI**K
         ENDIF
         IF(L.EQ.0.AND.PYF.EQ.0.D0) THEN
            PYFL=1.D0
         ELSE
            PYFL=PYF**L
         ENDIF

      XF = XF +
     &      J * A(I+1,J+1,K+1,L+1) *
     &      XII * PXFJ * YIK * PYFL
      ENDIF
             END DO
         END DO
          END DO
      END DO

      PXI=0.D0
      DO I=1,NORDNG
          DO J=0,NORDNG
         DO K=0,NORDNG
             DO L=0,NORDNG

      IF( I+J+K+L .LE. NORDNG ) THEN

         IF(I-1.EQ.0 .AND. XI.EQ. 0.D0)   THEN
            XII=1.D0
         ELSE
            XII=XI**(I-1)
         ENDIF
         IF(J.EQ.0 .AND. PXF.EQ. 0.D0)    THEN
            PXFJ=1.D0
         ELSE
            PXFJ=PXF**J
         ENDIF
         IF(K.EQ.0 .AND. YI.EQ.0.D0) THEN
            YIK=1.D0
         ELSE
            YIK=YI**K
         ENDIF
         IF(L.EQ.0.AND.PYF.EQ.0.D0) THEN
            PYFL=1.D0
         ELSE
            PYFL=PYF**L
         ENDIF

      PXI = PXI +
     &      I * A(I+1,J+1,K+1,L+1) *
     &      XII * PXFJ * YIK * PYFL
      ENDIF
             END DO
         END DO
          END DO
      END DO

      YF=0.D0
      DO I=0,NORDNG
          DO J=0,NORDNG
         DO K=0,NORDNG
             DO L=1,NORDNG

      IF( I+J+K+L .LE. NORDNG ) THEN

         IF(I.EQ.0 .AND. XI.EQ. 0.D0)  THEN
            XII=1.D0
         ELSE
            XII=XI**I
         ENDIF
         IF(J.EQ.0 .AND. PXF.EQ. 0.D0)    THEN
            PXFJ=1.D0
         ELSE
            PXFJ=PXF**J
         ENDIF
         IF(K.EQ.0 .AND. YI.EQ.0.D0) THEN
            YIK=1.D0
         ELSE
            YIK=YI**K
         ENDIF
         IF(L-1.EQ.0.AND.PYF.EQ.0.D0) THEN
            PYFL=1.D0
         ELSE
            PYFL=PYF**(L-1)
         ENDIF

      YF = YF +
     &      L * A(I+1,J+1,K+1,L+1) *
     &      XII * PXFJ * YIK * PYFL
      ENDIF
             END DO
         END DO
          END DO
      END DO

      PYI=0.D0
      DO I=0,NORDNG
          DO J=0,NORDNG
         DO K=1,NORDNG
             DO L=0,NORDNG
      IF( I+J+K+L .LE. NORDNG ) THEN

         IF(I.EQ.0 .AND. XI.EQ. 0.D0)  THEN
            XII=1.D0
         ELSE
            XII=XI**I
         ENDIF
         IF(J.EQ.0 .AND. PXF.EQ. 0.D0)    THEN
            PXFJ=1.D0
         ELSE
            PXFJ=PXF**J
         ENDIF
         IF(K-1.EQ.0 .AND. YI.EQ.0.D0) THEN
            YIK=1.D0
         ELSE
            YIK=YI**(K-1)
         ENDIF
         IF(L.EQ.0.AND.PYF.EQ.0.D0) THEN
            PYFL=1.D0
         ELSE
            PYFL=PYF**L
         ENDIF

      PYI = PYI +
     &      K * A(I+1,J+1,K+1,L+1) *
     &      XII * PXFJ * YIK * PYFL
      ENDIF
             END DO
         END DO
          END DO
      END DO

C--- TAKE CHANGE OF COORDINATES INTO ACCOUNT IN THE CALLING SEQUENCE

        CALL MYBFELD(ZI,YF,XF,BXI,BYI,BZI,AXI,AYI,AZI)

      AXI=AXI-AXF0
      AYI=AYI-AYF0
      AZI=AZI-AZF0

      DSQ1=1.D0-(PXI+AZI/BRHO)**2-(PYI+AYI/BRHO)**2

      IF (DSQ1.LE.0.) THEN

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     VECTOR POTENTIAL OR COEFFICIENTS CRAZY'
      WRITE(LUNGFO,*)'     1.D0-(PXI+AZI/BRHO)**2-(PYI+AYI/BRHO)**2 .LE.0'
      WRITE(LUNGFO,*)'     CHECK INPUT TO WAVE'
      WRITE(LUNGFO,*)

      WRITE(6,*)'     *** ERROR IN ERZFUN ***'
      WRITE(6,*)'     VECTOR POTENTIAL OR COEFFICIENTS CRAZY'
      WRITE(6,*)'     1.D0-(PXI+AZI/BRHO)**2-(PYI+AYI/BRHO)**2 .LE.0'
      WRITE(6,*)'     CHECK INPUT TO WAVE'
      STOP

      ENDIF

      DSQ1=1.D0/DSQRT(DSQ1)

C--- THE VECTOR POTENTIAL IS MORE OR LESS ARBITRARY BUT PXI,AZI AND XPI
C    MUST BE CONSISTENT

      XPI=(PXI+AZI/BRHO)*DSQ1
      YPI=(PYI+AYI/BRHO)*DSQ1

      RETURN
      END
