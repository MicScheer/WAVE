*CMZ :  0.00/06 16/02/2004  17.22.39  by  Michael Scheer
*CMZ :  0.00/05 23/12/2003  11.36.18  by  Michael Scheer
*CMZ :  0.00/04 19/12/2003  18.17.16  by  Michael Scheer
*CMZ :  0.00/02 12/12/2003  10.40.20  by  Michael Scheer
*CMZ :  0.00/01 08/12/2003  10.56.22  by  Michael Scheer
*-- Author :    Michael Scheer   04/12/2003
      subroutine bpen(imag,iplan,p1,p2,p3,vn)
c      subroutine bpen(vmag,p1,p2,p3,vn,idreh)

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
      double precision p1(3),p2(3),p3(3),vn(3),d2(3),d3(3),vnn
      integer imag,iplan
c      integer idreh,vmag(3)

c vmag is point inside magnet
c for the time being (23.12.2003), normal vector is calculated from direction
c of rotation


      d2(1)=p2(1)-p1(1)
      d2(2)=p2(2)-p1(2)
      d2(3)=p2(3)-p1(3)

      d3(1)=p3(1)-p1(1)
      d3(2)=p3(2)-p1(2)
      d3(3)=p3(3)-p1(3)

c vn= d2 x d3

      vn(1)=d2(2)*d3(3)-d2(3)*d3(2)
      vn(2)=d2(3)*d3(1)-d2(1)*d3(3)
      vn(3)=d2(1)*d3(2)-d2(2)*d3(1)

c      d2(1)=p1(1)-vmag(1)
c      d2(2)=p1(2)-vmag(2)
c      d2(3)=p1(3)-vmag(3)
c      vnn=vn(1)*d2(1)+vn(2)*d2(2)+vn(3)*d2(3)

        vnn=sqrt(vn(1)*vn(1)+vn(2)*vn(2)+vn(3)*vn(3))

c      if (vnn.lt.0.d0) THEN
c        vnn=-sqrt(vn(1)*vn(1)+vn(2)*vn(2)+vn(3)*vn(3))
c        idreh=-1
c      else
c        vnn=sqrt(vn(1)*vn(1)+vn(2)*vn(2)+vn(3)*vn(3))
c        idreh=1
c      endif

c      print *
c      print *,'idreh:', idreh
c      print *

      if (vnn.ne.0.d0) then
        vn(1)=vn(1)/vnn
        vn(2)=vn(2)/vnn
        vn(3)=vn(3)/vnn
      else
        print *,'***Error in bpen: Zero normal vector'
        print *,'magnet, plane:', imag,iplan
        stop
      endif

c      print *
c      print *,'vn:', vn(1),vn(2),vn(3)
c      print *

      return
      end
