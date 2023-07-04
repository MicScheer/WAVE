*CMZ :          28/11/2021  12.53.40  by  Michael Scheer
*-- Author :    Michael Scheer   23/08/2018
      subroutine util_lorentz_transformation(beta,gamma,boostmat,istat)

!+PATCH,//UTIL/FOR
!+DECK,util_lorentz_transformation.

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

      implicit none

      double precision :: beta(3),boostmat(4,4),bx,by,bz,gamma,b0,bn,
     &  g1,bxx,byy,bzz,bxy,bxz,byz,b2

      integer istat

      istat=0

      bn=norm2(beta)

      if (gamma.le.1.0d0) istat=-1

      if (bn.eq.0.0d0.or.gamma.le.1.0d0) then
        boostmat=0.0d0
        boostmat(1,1)=1.0d0
        boostmat(2,2)=1.0d0
        boostmat(3,3)=1.0d0
        boostmat(4,4)=1.0d0
        return
      endif

      b0=sqrt((1.0d0-1.0d0/gamma)*(1.0d0+1.0d0/gamma))
      b2=b0**2
      b0=b0/bn

      bx=beta(1)*b0
      by=beta(2)*b0
      bz=beta(3)*b0

      g1=gamma-1.0d0
      bxx=bx*bx/b2
      byy=by*by/b2
      bzz=bz*bz/b2
      bxy=bx*by/b2
      bxz=bx*bz/b2
      byz=by*bz/b2

      boostmat(1,1)=1.0d0+g1*bxx
      boostmat(1,2)=g1*bxy
      boostmat(1,3)=g1*bxz
      boostmat(1,4)=-gamma*bx

      boostmat(2,1)=boostmat(1,2)
      boostmat(2,2)=1.0d0+g1*byy
      boostmat(2,3)=g1*byz
      boostmat(2,4)=-gamma*by

      boostmat(3,1)=boostmat(1,3)
      boostmat(3,2)=boostmat(2,3)
      boostmat(3,3)=1.0d0+g1*bzz
      boostmat(3,4)=-gamma*bz

      boostmat(4,1)=boostmat(1,4)
      boostmat(4,2)=boostmat(2,4)
      boostmat(4,3)=boostmat(3,4)
      boostmat(4,4)=gamma

      return
      end
