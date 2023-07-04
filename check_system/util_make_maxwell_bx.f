*CMZ :          17/09/2018  11.42.31  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_make_maxwell_bx(nx,ny,nz,b,dx,dy,dz,bxoff,
     &  istatus)

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

c +PATCH,//UTIL/FOR
c +DECK,util_make_maxwell_bx

      implicit none

      include 'phyconparam.cmn'

      double precision, dimension (:,:,:), allocatable :: dbx
      double precision, dimension (:,:), allocatable :: dby,dbz
      double precision, dimension (:), allocatable :: by,bz,xs,ys,dbxs,
     &  ws1,ws2,ws3,ws4,cy,cz,zs,bys,bzs,bx,cx,cpy,cpz

      double precision, intent(inout) ::
     &  b(3,nx,ny,nz)
      double precision, intent(in) ::
     &  dx,dy,dz, bxoff

      integer, intent(in) :: nx,ny,nz
      integer, intent(out) :: istatus

      integer ix,iy,iz,nyz

      nyz=max(ny,nz)

      allocate(ys(ny),zs(nz),by(ny),bz(nz),cy(ny),cz(nz),
     &  ws1(nyz),ws2(nyz),ws3(nyz),ws4(nyz),dbz(ny,nz),dby(ny,nz),
     &  dbx(nx,ny,nz),dbxs(nx),bzs(nz),bys(ny),xs(nx),cpy(ny),cpz(nz))

      do iy=1,ny
        ys(iy)=iy*dy
      enddo
      do iz=1,nz
        zs(iz)=iz*dz
      enddo

      do ix=1,nx

        do iy=1,ny
          do iz=1,nz
            bzs(iz)=b(3,ix,iy,iz)
          enddo
          call util_spline_coef_deriv(zs,bzs,nz,9999.0d0,9999.0d0,
     &      cpz,cz,ws1,ws2,ws3,ws4)
          dbz(iy,1:nz)=cpz(1:nz)
        enddo

        do iz=1,nz
          bys(1:ny)=b(2,ix,1:ny,iz)
          call util_spline_coef_deriv(ys,bys,ny,9999.0d0,9999.0d0,
     &      cpy,cy,ws1,ws2,ws3,ws4)
          dby(1:ny,iz)=cpy(1:ny)
        enddo

        do iz=1,nz
          do iy=1,ny
            dbx(ix,iy,iz)=-dby(iy,iz)-dbz(iy,iz)
          enddo
        enddo

        xs(ix)=ix*dx

      enddo !nx

      deallocate(ws1,ws2,ws3,ws4)
      allocate(bx(nx),ws1(nx),ws2(nx),ws3(nx),ws4(nx),cx(nx))

      do iz=1,nz
          do iy=1,ny
            dbxs(1:nx)=dbx(1:nx,iy,iz)
            call util_spline_running_integral(xs,dbxs,nx,bx,cx,
     &        ws1,ws2,ws3,ws4)
            b(1,1:nx,iy,iz)=bx(1:nx)+bxoff
          enddo
        enddo

      istatus=0

c9999  continue

      deallocate(ys,zs,by,bz,cy,cz,
     &  ws1,ws2,ws3,ws4,dbz,dby,
     &  dbx,dbxs,bx,bzs,bys,xs,cpy,cpz)


      return
      end
