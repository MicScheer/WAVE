*CMZ :  4.00/15 29/04/2022  11.54.49  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.15/00 01/05/2000  13.23.51  by  Michael Scheer
*CMZ : 00.01/07 15/03/95  15.04.35  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_MAX_PARABEL(NDIM,X,F,XMX,FMX,WSX,WSF,JFAIL)
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

C--- TO FIND MAXIMUM OF ARRAY FUNCTION F(X)

C     INPUT : F   ARRAY OF FUNCTION
C             X   ARRAY OF ARGUMENTS
C             WSX WORKINGSPACE
C             WSF WORKINGSPACE

C     OUTPUT:  XMX ARGUMENT WHERE FUNCTION REACHES EXTREMUM
C              FMX MAXIMUM OF FUNCTION
C              JFAIL FLAG: =0, IF OK, =1 ELSE

      IMPLICIT NONE

      INTEGER NDIM,I,IFAIL,JFAIL,IMAX,IFOUND
      DOUBLE PRECISION X(NDIM),F(NDIM),XMX,FMX,WSX(NDIM),WSF(NDIM)
      DOUBLE PRECISION XDUM(3),FDUM(3),A(3),FP(3),FMAX,XMAX

      IFOUND=0
      IMAX=0
      FMAX=-1.0D30

      if (ndim.le.0) then
        jfail=-1
        print*, "*** Error in util_max_parabel: Ndim .le. 0"
        return
      endif

      DO I=1,NDIM
        WSX(I)=X(I)
        WSF(I)=F(I)
        IF (WSF(I).NE.0.0D0) IFOUND=1
        IF (WSF(I).GT.FMAX) THEN
          FMAX=WSF(I)
          XMAX=WSX(I)
          IMAX=I
        ENDIF
      ENDDO

      if (imax.eq.0) then
        jfail=-2
        print*, "*** Error in util_max_parabel: IMAX = 0"
        return
      endif
C      CALL UTIL_SORT_FUNC(NDIM,WSF,WSX)
C      DO I=1,3
C         XDUM(I)=WSX(NDIM-I+1)
C         FDUM(I)=WSF(NDIM-I+1)
C      ENDDO

      IF (IFOUND.EQ.0) THEN
        XMX=0.D0
        FMX=0.D0
        JFAIL=0
        RETURN
      ENDIF

      IF (IMAX.EQ.1.OR.IMAX.EQ.NDIM) THEN
          JFAIL=1
          FMX=FMAX
          XMX=XMAX
          RETURN
      ENDIF

        XDUM(1)=WSX(IMAX-1)
        FDUM(1)=WSF(IMAX-1)
        XDUM(2)=WSX(IMAX)
        FDUM(2)=WSF(IMAX)
        XDUM(3)=WSX(IMAX+1)
        FDUM(3)=WSF(IMAX+1)

      CALL UTIL_PARABEL(XDUM,FDUM,A,FP,XMX,FMX,IFAIL)

      IF (IFAIL.NE.0.OR.FMX.LT.FMAX) THEN
          JFAIL=1
          FMX=FMAX
          XMX=XMAX
          RETURN
      ENDIF

      JFAIL=0

      RETURN
      END
