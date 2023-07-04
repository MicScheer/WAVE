*CMZ :  2.52/14 20/12/2004  09.48.45  by  Michael Scheer
*-- Author :    Michael Scheer   17/12/2004
      subroutine corner_mat(eps,rho,hmat,vmat,b)
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

      double precision eps,rho,vmat(2,2),hmat(2,2),ts,cs,b

      cs=cos(eps)

      if (cs.eq.0.0d0) then
        print*,'*** Error in corner_mat: Cosine of corner angle is zero'
        stop '*** Program WAVE aborted ***'
      endif

      ts=tan(eps)

      if (rho.ne.0.0d0) then

        hmat(1,1)=1.0d0
        hmat(1,2)=0.0d0
        hmat(2,1)=ts/rho
        hmat(2,2)=1.0d0

        vmat(1,1)=1.0d0
        vmat(1,2)=0.0d0
        vmat(2,1)=(b/(6.0d0*rho*cs)-ts)/rho
        vmat(2,2)=1.0d0

      else

        hmat(1,1)=1.0d0
        hmat(1,2)=rho !we assume rho to be length of a drift, here
        hmat(2,1)=0.0d0
        hmat(2,2)=1.0d0

        vmat=hmat

      endif

      return
      end
