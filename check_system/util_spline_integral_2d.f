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

      double precision x(nx),y(ny),f(nx,ny),result
      integer istat,nx,ny,ix,iy,kstat

      double precision, allocatable :: fb(:),fb2(:),coef(:),
     &  w1(:),w2(:),w3(:),w4(:)

      allocate(fb(max(nx,ny)))
      allocate(fb2(max(nx,ny)))
      allocate(coef(max(nx,ny)))
      allocate(w1(max(nx,ny)))
      allocate(w2(max(nx,ny)))
      allocate(w3(max(nx,ny)))
      allocate(w4(max(nx,ny)))

      kstat=0

      if (ny.gt.nx) then
        do ix=1,nx
          do iy=1,ny
            fb(iy)=f(ix,iy)
            call util_spline_integral_stat(y(iy),fb,ny,fb2(ix)
     &        ,coef,w1,w2,w3,w4,istat)
            kstat=kstat+istat
          enddo
          call util_spline_integral_stat(x(ix),fb2,nx,result
     &      ,coef,w1,w2,w3,w4,istat)
          kstat=kstat+istat
        enddo
      else !nx.gt.ny?
        do iy=1,ny
          do ix=1,nx
            fb(ix)=f(ix,iy)
            call util_spline_integral_stat(x(ix),fb,nx,fb2(iy)
     &        ,coef,w1,w2,w3,w4,istat)
            kstat=kstat+istat
          enddo
          call util_spline_integral_stat(y(iy),fb2,ny,result
     &      ,coef,w1,w2,w3,w4,istat)
          kstat=kstat+istat
        enddo
      endif !nx.gt.ny

      deallocate(fb)
      deallocate(fb2)
      deallocate(coef)
      deallocate(w1)
      deallocate(w2)
      deallocate(w3)
      deallocate(w4)

      return
      end
