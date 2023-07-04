*CMZ :  4.00/11 07/05/2021  07.34.04  by  Michael Scheer
*CMZ :  4.00/07 09/01/2020  13.45.07  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.29.14  by  Michael Scheer
*CMZ :  2.14/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  12.40.42  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.09.58  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  16.56.14  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.51.02  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.55  by  Michael Scheer
*-- Author : Chaoen Wang
C...............................................................................
      Subroutine FFT(Z1,M,ISign)
C...............................................................................
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

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      COMPLEX Z1(NDOBSVZP+NDOBSVYP) !SR USMCON2

CMSH  Parameter (MaxPt0=616)
CMSH  CompLex Z1(MaxPt0)

      CompLex ZU,ZW,ZT

      INTEGER K,NM1,NV2,I,J,N,IP,LE,LE1,L,M,ISIGN

      REAL PI

      DATA PI/3.141592653589793D0/

      print*,"FFT IS OBSOLETE, SEE //WAVE/USEM IN WAVE.CMZ"

      Return
      End
