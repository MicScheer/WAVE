*CMZ :          03/08/2025  10.05.09  by  Michael Scheer
*CMZ :  4.01/03 17/05/2023  11.24.58  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  11.49.33  by  Michael Scheer
*CMZ : 00.00/16 21/11/2014  14.53.59  by  Michael Scheer
*-- Author :    Michael Scheer   21/11/2014
      subroutine util_spline_integral_2d(nx,ny,x,y,f,result,istat)
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

      double precision x(nx),y(ny),f(nx,ny),result,coef(4,4,nx,ny),a(4,4),dx,dy
      integer nx,ny,ix,iy,istat,k,l

      call util_coef_spline_2d(nx,ny,f,coef,istat)
      if (istat.ne.0) return

      dx=(x(nx)-x(1))/(nx-1)
      dy=(y(ny)-y(1))/(ny-1)

      result=0.0d0

      do iy=1,ny-1
        do ix=1,nx-1

          a(1:4,1:4)=coef(1:4,1:4,ix,iy)
          do k=1,4
            do l=1,4
              result=result+a(k,l)/k/l
            enddo
          enddo

        enddo
      enddo

c      a(1:4,1:4)=coef(1:4,1:4,1,1)
c      do k=1,4
c        do l=1,4
c          result=result-a(k,l)/k/l
c        enddo
c      enddo

      result=result*dx*dy

      end
