*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  16.22.05  by  Michael Scheer
*CMZ :  1.00/00 09/06/97  12.17.43  by  Michael Scheer
*-- Author : Michael Scheer   09/06/97

      SUBROUTINE BLINE(X,Y,Z,BX,BY,BZ)
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

      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,CURRU,CENXU,CENYU,CURRD,CENXD,CENYD
     &,RX1,RX2,RY1,RY2,R31,R32
     &,RX3,RX4,RY3,RY4,R33,R34
     &,BXU,BXD,BYU,BYD,BZU,BZD

*KEEP,uservar.
      include 'uservar.cmn'
*KEND.

      Z=Z

      CURRU=USER(1)
      CENXU=USER(2)
      CENYU=USER(3)
      CURRD=USER(4)
      CENXD=USER(5)
      CENYD=USER(6)

      RX1=X-CENXU
      RY1=Y-CENYU
      RX2=X+CENXU
      RY2=Y-CENYU

      RX3=X-CENXD
      RY3=Y-CENYD
      RX4=X+CENXD
      RY4=Y-CENYD

      R31=DABS(RX1*RX1+RY1*RY1)
      R32=DABS(RX2*RX2+RY2*RY2)
      R33=DABS(RX3*RX3+RY3*RY3)
      R34=DABS(RX4*RX4+RY4*RY4)

      BXU= CURRU*(RY1/R31-RY2/R32)
      BYU=-CURRU*(RX1/R31-RX2/R32)
      BZU=0.D0
      BXD= CURRD*(RY3/R33-RY4/R34)
      BYD=-CURRD*(RX3/R33-RX4/R34)
      BZD=0.D0

      BX=BXU+BXD
      BY=BYU+BYD
      BZ=BZU+BZD

      END
