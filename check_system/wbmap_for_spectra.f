*CMZ :  4.00/15 09/03/2022  15.57.19  by  Michael Scheer
*-- Author :    Michael Scheer   28/09/95
      subroutine wbmap_for_spectra
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

c--- to write 3d field map for program spectra

      implicit none

      external dcosd,dsind
      double precision dcosd,dsind

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bmap.
      include 'bmap.cmn'
*KEND.

      double precision :: dx=0.0d0,dy=0.0d0,dz=0.0d,
     &  bx,by,bz,ax,ay,az,x,y,z,r,phi

      integer ix,iy,iz,lunmap

      open(newunit=lunmap,file=filebmap,status='new')

      if (xmapmn.eq.9999.) xmapmn=xstart
      if (xmapmx.eq.9999.) xmapmx=xstop

      if (xmapmx.lt.xmapmn) then
        write(lungfo,*)
        write(lungfo,*)'*** Error in wbmap_to_spectra xmapmx.lt.xmapmn'
        write(lungfo,*)'Check namelist wbmap in wave.in'
        write(lungfo,*)
        write(6,*)
        write(6,*)'*** Error in wbmap_to_spectra xmapmx.lt.xmapmn'
        write(6,*)'Check namelist wbmap in wave.in'
        write(6,*)
        stop '*** Program WAVE aborted ***'
      endif

      if (ymapmx.lt.ymapmn) then
        write(lungfo,*)
        write(lungfo,*)'*** error in wbmap_to_spectra ymapmx.lt.ymapmn'
        write(lungfo,*)'check namelist wbmap in wave.in'
        write(lungfo,*)
        write(6,*)
        write(6,*)'*** error in wbmap_to_spectra ymapmx.lt.ymapmn'
        write(6,*)'check namelist wbmap in wave.in'
        write(6,*)
        stop '*** Program WAVE aborted ***'
      endif

      if (zmapmx.lt.zmapmn) then
        write(lungfo,*)
        write(lungfo,*)'*** error in wbmap_to_spectra zmapmx.lt.zmapmn'
        write(lungfo,*)'check namelist wbmap in wave.in'
        write(lungfo,*)
        write(6,*)
        write(6,*)'*** error in wbmap_to_spectra zmapmx.lt.zmapmn'
        write(6,*)'check namelist wbmap in wave.in'
        write(6,*)
        stop '*** Program WAVE aborted ***'
      endif

      if (nmapx.eq.-9999) nmapx=nint((xstop-xstart)*myinum)+1

      if (nmapx.gt.1) dx=(xmapmx-xmapmn)/(nmapx-1)
      if (nmapy.gt.1) dy=(ymapmx-ymapmn)/(nmapy-1)
      if (nmapz.gt.1) dz=(zmapmx-zmapmn)/(nmapz-1)

      write(lunmap,*) sngl(dz*1000.0d0),sngl(dy*1000.0d0),sngl(dx*1000.0d0),
     &  nmapz,nmapy,nmapx

      do iz=1,nmapz
        do iy=1,nmapy
          do ix=1,nmapx
            x=xmapmn+(ix-1)*dx
            y=ymapmn+(iy-1)*dy
            z=-(zmapmn+(iz-1)*dz)
            call mybfeld(x,y,z,bx,by,bz,ax,ay,az)
            write(lunmap,*)sngl(bx),sngl(by),sngl(-bz)
            enddo
          enddo
        enddo

      flush (lunmap)
      close (lunmap)

      write(lungfo,*)
      write(lungfo,*)'     Subroutine wbmap_to_spectra: field map written'
      write(lungfo,*)'       File:        ',filebmap
      write(lungfo,*)'       nx,xmin,xmax:',nmapx,xmapmn,xmapmx
      write(lungfo,*)'       ny,ymin,ymax:',nmapy,ymapmn,ymapmx
      write(lungfo,*)'       nz,zmin,zmax:',nmapz,zmapmn,zmapmx
      write(lungfo,*)

      return
      end
