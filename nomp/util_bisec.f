*CMZ :  2.63/03 15/05/2009  11.26.16  by  Michael Scheer
*CMZ :  2.63/02 01/02/2008  16.04.42  by  Michael Scheer
*CMZ : 00.00/02 25/08/2006  15.27.06  by  Michael Scheer
*CMZ : 00.00/01 23/02/96  14.56.50  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.54  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_bisec(n,xa,x,i)
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

c returns i with:
c x.ge.x(i) .and. x.lt.x(i+1), if x(1).lt.x(n)
c or
c x.le.x(i) .and. x.gt.x(i+1), if x(n).lt.x(1)

      IMPLICIT NONE

      REAL*8 XA(*),x
      integer n,i,klo,khi,k

      i=-1

      if (n.lt.2) then
        return
      else if (x.eq.xa(1)) then
        i=1
        return
      else if (x.eq.xa(n)) then
        i=n-1
        return
      endif

      if (xa(1).lt.xa(n)) then

        if (x.lt.xa(1).or.x.gt.xa(n)) return

        klo=1
        KHI=N

1       IF (KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XA(K).GT.X)THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
          GOTO 1
        ENDIF

      else if (xa(n).lt.xa(1)) then

        if (x.lt.xa(n).or.x.gt.xa(1)) return

        klo=1
        Khi=n

2       IF (KHI-KLO.GT.1) THEN

          K=(KHI+KLO)/2

          IF(XA(K).LT.X)THEN
            Khi=K
          ELSE
            Klo=K
          ENDIF

          GOTO 2

        ENDIF

      endif !xa(1).lt.xa(n)

      i=klo

      RETURN
      END
