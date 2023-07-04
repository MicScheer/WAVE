*CMZ :  1.00/00 29/07/97  10.17.24  by  Michael Scheer
*CMZ : 00.02/00 15/11/96  11.30.41  by  Michael Scheer
*CMZ : 00.01/09 25/10/95  17.37.17  by  Michael Scheer
*-- Author :    Michael Scheer   29/09/95

      SUBROUTINE BMASHB(X,Y,Z,BX,BY,BZ)
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

C--- ALL INDICES ACCORDING TO FORTRAN, BUT LORD3DG AND MORD3DG REFER TO MATH. INDICES

      IMPLICIT NONE

*KEEP,bmash.
      include 'bmash.cmn'
*KEND.


      INTEGER IX,IY,IZ,IND

      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,XPOW1,YPOW1,ZPOW1

      DOUBLE PRECISION XPOW(MORDP+1),YPOW(MORDP+1),ZPOW(MORDP+1)


      XPOW(1)=1.D0
      DO IX=2,MORDP+1
          XPOW(IX)=XPOW(IX-1)*X
      ENDDO

      YPOW(1)=1.D0
      DO IY=2,MORDP+1
          YPOW(IY)=YPOW(IY-1)*Y
      ENDDO

      ZPOW(1)=1.D0
      DO IZ=2,MORDP+1
          ZPOW(IZ)=ZPOW(IZ-1)*Z
      ENDDO

      BX=0.D0
      BY=0.D0
      BZ=0.D0

      DO IND=1,MTOTP

             IX=ICINDG(1,IND)
             IY=ICINDG(2,IND)
             IZ=ICINDG(3,IND)

             IF (IX.GT.1) THEN
            XPOW1=XPOW(IX-1)
             ELSE
            XPOW1=1.D0
             ENDIF
             BX=BX-(IX-1)*CINDG(IND)*XPOW1*YPOW(IY)*ZPOW(IZ)

             IF (IY.GT.1) THEN
            YPOW1=YPOW(IY-1)
             ELSE
            YPOW1=1.D0
             ENDIF

             BY=BY-(IY-1)*CINDG(IND)*XPOW(IX)*YPOW1*ZPOW(IZ)

             IF (IZ.GT.1) THEN
            ZPOW1=ZPOW(IZ-1)
             ELSE
            ZPOW1=1.D0
             ENDIF

             BZ=BZ-(IZ-1)*CINDG(IND)*XPOW(IX)*YPOW(IY)*ZPOW1

      ENDDO

      RETURN
      END
