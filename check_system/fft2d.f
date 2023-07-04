*CMZ :  4.00/07 09/01/2020  13.44.20  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/08 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.14/02 19/04/2000  17.02.45  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  12.15.08  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.09.58  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  16.56.43  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.51.06  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.55  by  Michael Scheer
*-- Author : Chaoen Wang
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
C..........................................................................
        Subroutine FFT2d(NoMx,NoMy,HXLgh,HYLgh,ISign)
C..........................................................................

*KEEP,usemf90u.
      include 'usemf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      COMPLEX Z1(NDOBSVZP+NDOBSVYP) !SR USMCON2

CMSH  Parameter (MaxPt0=616,MaxPtx=616,MaxPty=553)

CMSH        CompLex Z1(MaxPt0),Z2(MaxPtx,MaxPty)

        CompLex ZI,ZJ

      INTEGER NODELX,NODELY,NOMX,NOMY,ISIGN,JYTH,JXTH,NYTH,IXTH,MXTH
      REAL HKyLgh,DeLKy,HKxLgh,DeLy,DeLkx,DeLx,HXLgh,HYLgh,PI

      DATA PI/3.141592653589793D0/

      PRINT*,"FFT2D IS OBSOLETE, SEE //WAVE/USEM IN WAVE.CMZ"
        Return
      End
