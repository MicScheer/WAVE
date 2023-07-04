*CMZ :  0.00/06 14/01/2004  15.36.48  by  Michael Scheer
*CMZ :  0.00/05 23/12/2003  14.17.24  by  Michael Scheer
*CMZ :  0.00/04 23/12/2003  10.12.22  by  Michael Scheer
*CMZ :  0.00/03 15/12/2003  14.19.52  by  Michael Scheer
*CMZ :  0.00/02 12/12/2003  10.32.59  by  Michael Scheer
*CMZ :  0.00/01 08/12/2003  11.19.50  by  Michael Scheer
*-- Author :    Michael Scheer   04/12/2003
      subroutine bpet(vnz,t,t1)

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
      double precision vnz(3),t(3,3),t1(3,3)
      double precision det,eps,vz1,vxyz1,vyyz1

      integer i,j
      integer ifail

      data eps/1.0d-10/

      if (abs(vnz(3)+1.0d0).lt.eps) then

        t(1,1)=-1.d0
        t(1,2)=0.d0
        t(1,3)=0.d0

        t(2,1)=0.d0
        t(2,2)=1.d0
        t(2,3)=0.d0

        t(3,1)=0.d0
        t(3,2)=0.d0
        t(3,3)=-1.d0

      else

        vz1=vnz(3)+1.d0
        vxyz1=-vnz(1)*vnz(2)/vz1
        vyyz1=vnz(2)*vnz(2)/vz1

        t(1,1)=vnz(3)+vyyz1
        t(1,2)=vxyz1
        t(1,3)=-vnz(1)

        t(2,1)=vxyz1
        t(2,2)=1.d0-vyyz1
        t(2,3)=-vnz(2)

        t(3,1)=vnz(1)
        t(3,2)=vnz(2)
        t(3,3)=vnz(3)

      endif

      do i=1,3
        do j=1,3
          t1(i,j)=t(j,i)
        enddo
      enddo

      call util_determinante(3,t,det,ifail)

      if (ifail.ne.0) then
        print *,'*** Error in bpet: Bad determinante of matrix'
        stop
      endif

      if (abs(det-1.d0).gt.1.d-10) then
        print *,'*** Error in bpet: Bad determinante of matrix'
        print *
        print *,'Matrix t'
        do i=1,3
          print *,t(i,1),t(i,2),t(i,3)
        enddo
        print *
        print *,'det: ',det
        print *
        stop
      endif


      return
      end
