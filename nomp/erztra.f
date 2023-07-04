*CMZ :  3.00/00 11/03/2013  10.37.07  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.50.22  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.54  by  Michael Scheer
*-- Author : Michael Scheer
C****************************************************************
      SUBROUTINE ERZTRA (XI,PXF,YI,PYF,   XF,PXI,YF,PYI)
C****************************************************************
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
C
C FUNCTIONAL DESCRIPTION:
C
C     CALCULATES CANONICAL COORDINATES AND MOMENTA ACCORDING TO
C     COEFFICIENTS ADUM (I,J,K,L) OF GENERATING FUNCTION
C
C     INPUT: XI,PXF,YI,PYF, ADUM VIA COMMOM/AKO/
C     OUTPUT:XF,PXI,YF,PYI
C
      IMPLICIT NONE

      INTEGER I,J,K,L,NORD,NKOEF

*KEEP,genfun.
      include 'genfun.cmn'
*KEND.

      DOUBLE PRECISION  XI,PXF,YI,PYF,XF,PXI,YF,PYI,XII,PXFJ,YIK,PYFL

      NORD=NORDNG-1

      XF=0.D0
      DO I=0,NORD
          DO J=1,NORD
         DO K=0,NORD
             DO L=0,NORD

      IF( I+J+K+L .LE. NORD ) THEN

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
     &      J * ADUM(I+1,J+1,K+1,L+1) *
     &      XII * PXFJ * YIK * PYFL
      ENDIF
             END DO
         END DO
          END DO
      END DO

      PXI=0.D0
      DO I=1,NORD
          DO J=0,NORD
         DO K=0,NORD
             DO L=0,NORD

      IF( I+J+K+L .LE. NORD ) THEN

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
     &      I * ADUM(I+1,J+1,K+1,L+1) *
     &      XII * PXFJ * YIK * PYFL
      ENDIF
             END DO
         END DO
          END DO
      END DO

      YF=0.D0
      DO I=0,NORD
          DO J=0,NORD
         DO K=0,NORD
             DO L=1,NORD

      IF( I+J+K+L .LE. NORD ) THEN

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
     &      L * ADUM(I+1,J+1,K+1,L+1) *
     &      XII * PXFJ * YIK * PYFL
      ENDIF
             END DO
         END DO
          END DO
      END DO

      PYI=0.D0
      DO I=0,NORD
          DO J=0,NORD
         DO K=1,NORD
             DO L=0,NORD
      IF( I+J+K+L .LE. NORD ) THEN

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
     &      K * ADUM(I+1,J+1,K+1,L+1) *
     &      XII * PXFJ * YIK * PYFL
      ENDIF
             END DO
         END DO
          END DO
      END DO

      RETURN
      END
