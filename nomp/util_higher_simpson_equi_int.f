*CMZ :  2.63/05 13/08/2009  15.38.27  by  Michael Scheer
*-- Author :    Michael Scheer   13/08/2009
      subroutine util_higher_simpson_equidist_integral(n,x,f,sum,istat)
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

      double precision x(*),f(*),sum,sumeve,sumodd,dx
      integer n,i,istat

      istat=-1

      if (n.lt.2) return

      sum=0.0d0
      sumeve=0.0d0
      sumodd=0.0d0

      dx=(x(n)-x(1))/(n-1)

      if (n.ge.8) then
        do i=5,n-4
          sumeve=sumeve+f(i)
        enddo
        sum=(
     &    (17.0d0*f(1)+59.0d0*f(2)+43.0d0*f(3)+49.0d0*f(4))/48.0d0
     &   +sumeve
     &    +(17.0d0*f(n)+59.0d0*f(n-1)+43.0d0*f(n-2)+49.0d0*f(n-3))/48.0d0
     &    )*dx
      else if (n.eq.7) then
        sum=(f(1)+4.0d0*(f(2)+f(4)+f(6))+2.0d0*(f(3)+f(5))+f(7))*dx/3.0d0
      else if (n.eq.6) then
        sum=(f(3)+f(4)+(5.0d0*f(1)+13.0d0*f(2)
     &    +13.0d0*f(5)+5.0d0*f(6))/12.0d0)*dx
      else if (n.eq.5) then
        sum=(f(1)+4.0d0*(f(2)+f(4))+2.0d0*f(3)+f(5))*dx/3.0d0
      else if (n.eq.4) then
        sum=(5.0d0*f(1)+13.0d0*f(2)
     &    +13.0d0*f(3)+5.0d0*f(4))*dx/12.0d0
      else if (n.eq.3) then
        sum=(f(1)+4.0d0*f(2)+f(3))*dx/3.0d0
      else if (n.eq.2) then
        sum=(f(1)+f(2))*dx/2.0d0
      endif

      istat=0

      return
      end
