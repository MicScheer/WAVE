*CMZ :  2.66/09 22/03/2010  15.00.25  by  Michael Scheer
*CMZ :  2.66/07 10/12/2009  10.40.02  by  Michael Scheer
*CMZ : 00.00/01 20/06/95  10.09.17  by  Michael Scheer
*-- Author :    Michael Scheer   26/01/95
      SUBROUTINE util_min_parabel(NDIM,X,F,XMN,FMN,WSX,WSF,JFAIL)
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

C--- TO FIND MINIMUM OF ARRAY FUNCTION F(X)

C     INPUT : F(NDIM)   ARRAY OF FUNCTION
C             X(NDIM)   ARRAY OF ARGUMENTS
C             WSX(NDIM) WORKINGSPACE
C             WSF(NDIM) WORKINGSPACE

C     OUTPUT:  XMN ARGUMENT WHERE FUNCTION REACHES EXTREMUM
C              FMN MINIMUM OF FUNCTION
C              JFAIL FLAG: =0, IF OK, =1 ELSE

      IMPLICIT NONE

      INTEGER NDIM,I,IFAIL,JFAIL

      REAL*8 X(NDIM),F(NDIM),XMN,FMN,WSX(NDIM),WSF(NDIM)
      REAL*8 XDUM(3),FDUM(3),A(3),FP(3),FMIN,XMIN

      FMIN=1.0D30
      DO I=1,NDIM
          WSX(I)=X(I)
          WSF(I)=F(I)
          IF (WSF(I).LT.FMIN) THEN
              FMIN=WSF(I)
              XMIN=WSX(I)
          ENDIF
      ENDDO

      CALL UTIL_SORT_FUNC(NDIM,WSF,WSX)

      DO I=1,3
         XDUM(I)=WSX(I)
         FDUM(I)=WSF(I)
      ENDDO

      CALL UTIL_PARABEL(XDUM,FDUM,A,FP,XMN,FMN,IFAIL)

      IF (IFAIL.NE.0.OR.FMN.LT.FMIN) THEN
          JFAIL=1
c          WRITE(6,*)
c          WRITE(6,*)'*** WARNING IN UTIL_MIN_PARABEL: SEARCH FAILED ***'
c          WRITE(6,*)'*** MINIMUM OF ARRAY TAKEN ***'
c          WRITE(6,*)
          FMN=FMIN
          XMN=XMIN
          RETURN
      ENDIF

      JFAIL=0

      RETURN
      END
