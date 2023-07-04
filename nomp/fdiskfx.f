*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.14/02 19/04/2000  17.21.02  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  15.49.52  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.50.57  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.55  by  Michael Scheer
*-- Author : Michael Scheer
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
C...............................................................................
      Function FDiskFx(Rkx,Rky,Radius)
C...............................................................................

      IMPLICIT NONE

      REAL Rkx,Rky,Radius,FDISKFX,PI,BESJ1
      ReaL Dummy
      DATA PI/3.141592653589793D0/

CMSH  FDiskFx=1.0 !normierte
      FDiskFx=(Pi*Radius*Radius) !MSH
C     FDiskFx=(Pi*Radius*Radius)/(Pi*Radius*Radius) !normierte
C     FDiskFx=(Pi*Radius*Radius)
      If ( Rkx*Rkx+Rky*Rky.ne.0 ) Then
           Dummy=Radius*Sqrt(Rkx*Rkx+Rky*Rky)

CMSH USE CERN ROUTINE BESJ1 (C312)
           FDiskFx=FDiskFx*2.0*BESJ1(Dummy)/Dummy
CMSH       FDiskFx=FDiskFx*2.0*BesselFx1(Dummy)/Dummy

      Endif
        Return
        End
