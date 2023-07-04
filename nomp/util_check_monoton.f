*CMZ :  2.52/09 26/10/2004  11.55.46  by  Michael Scheer
*-- Author :    Michael Scheer   26/10/2004
      SUBROUTINE UTIL_CHECK_MONOTON(N,X,IMONO)
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

C     Returns IMONO: 2 for strictly raising,
C                    1 for raising,
C                    0 not monoton or const.
C                   -1 for falling,
C                   -2 for strictly falling,


      IMPLICIT NONE

      INTEGER I,N,IMONO,IMONO0,IUP
      DOUBLE PRECISION X(N)

      IUP=0
      DO I=1,N
        IF (X(1).LT.X(I)) THEN
          IUP=1
          GOTO 1
        ELSE IF (X(1).GT.X(I)) THEN
          IUP=-1
          GOTO 1
        ENDIF
      ENDDO

1     IMONO0=0
      IMONO=0

      IF (IUP.EQ.1) THEN

        DO I=1,N-1
          IF (X(I).EQ.X(I+1)) IMONO0=1
          IF (X(I).GT.X(I+1)) RETURN
        ENDDO

        IF (IMONO0.EQ.0) THEN
          IMONO=2
        ELSE
          IMONO=1
        ENDIF

      ELSE IF (IUP.EQ.-1) THEN

        DO I=1,N-1
          IF (X(I).EQ.X(I+1)) IMONO0=1
          IF (X(I).LT.X(I+1)) RETURN
        ENDDO

        IF (IMONO0.EQ.0) THEN
          IMONO=-2
        ELSE
          IMONO=-1
        ENDIF

      ELSE

        IMONO=0

      ENDIF

      RETURN
      END
