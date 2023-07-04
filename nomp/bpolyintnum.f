*CMZ :  3.05/00 25/04/2018  13.09.51  by  Michael Scheer
*CMZ :  2.52/05 16/08/2004  13.37.21  by  Michael Scheer
*CMZ :  1.02/00 28/07/2004  15.02.40  by  Michael Scheer
*-- Author :    Michael Scheer   27/07/2004
      subroutine bpolyintnum(kmag,xint,yint,zint,
     &  vxint,vyint,vzint,
     &  bxint,byint,bzint)
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

c calculates numerically first magnetic integral for line defined by point
c      (xint,yint,zint) and vector (vxint,vyint,vzint) in the range
c      ranginti->ranginte with nstepint steps

c xint,yint,zint given in meter
c bxint,byint,bzint given in Tm

*KEEP,bwpolyederf90u.
      include 'bwpolyederf90u.cmn'
*KEND.

      implicit none

      double precision xint,yint,zint,vxint,vyint,vzint,vnx,vny,vnz,vn,
     &  bxint,byint,bzint,xx,yy,zz,ds,ds2,bxi(4),byi(4),bzi(4),tiny2,x,y,z,
     &  bx,by,bz

      integer ical,i
      integer kmag

      data ical/0/

      if (ical.eq.0) then
        tiny2=tiny*tiny
        ical=1
      endif !ical

      vn=sqrt(vxint*vxint+vyint*vyint+vzint*vzint)

      vnx=vxint/vn
      vny=vyint/vn
      vnz=vzint/vn

      bxint=0.d0
      byint=0.d0
      bzint=0.d0

      bxi=0.0d0
      byi=0.0d0
      bzi=0.0d0

      if (nstepint.le.0.or.bpebc(7,kmag).eq.0.0d0) return

      nstepint=max(nstepint,4)/2*2+1
      ds=(ranginte-ranginti)/(nstepint-1)*0.001d0
      ds2=2.0d0*ds

      if (vnx.lt.tiny2.and.vny.lt.tiny2) then

        x=xint
        y=yint

        z=ranginti*0.001d0
        call bpolyeder1(kmag,x,y,z,bx,by,bz)
        bxi(1)=bx/3.0d0
        byi(1)=by/3.0d0
        bzi(1)=bz/3.0d0

        z=ranginte*0.001d0
        call bpolyeder1(kmag,x,y,z,bx,by,bz)
        bxi(4)=bx/3.0d0
        byi(4)=by/3.0d0
        bzi(4)=bz/3.0d0

        z=ranginti*0.001d0-ds
        do i=2,nstepint-1,2
          z=z+ds2
          call bpolyeder1(kmag,x,y,z,bx,by,bz)
          bxi(2)=bxi(2)+bx
          byi(2)=byi(2)+by
          bzi(2)=bzi(2)+bz
        enddo

        z=ranginti*0.001d0
        do i=3,nstepint-2,2
          z=z+ds2
          call bpolyeder1(kmag,x,y,z,bx,by,bz)
          bxi(3)=bxi(3)+bx
          byi(3)=byi(3)+by
          bzi(3)=bzi(3)+bz
        enddo

        bxint=
     &    (bxi(1)/3.0d0+4.0d0/3.0d0*bxi(2)+2.0d0/3.0d0*bxi(3)+bxi(4)/3.0d0)*ds
        byint=
     &    (byi(1)/3.0d0+4.0d0/3.0d0*byi(2)+2.0d0/3.0d0*byi(3)+bxi(4)/3.0d0)*ds
        bzint=
     &    (bzi(1)/3.0d0+4.0d0/3.0d0*bzi(2)+2.0d0/3.0d0*bzi(3)+bxi(4)/3.0d0)*ds

      else if (vnx.lt.tiny2.and.vnz.lt.tiny2) then

        x=xint
        z=zint

        y=ranginti*0.001d0
        call bpolyeder1(kmag,x,y,z,bx,by,bz)
        bxi(1)=bx/3.0d0
        byi(1)=by/3.0d0
        bzi(1)=bz/3.0d0

        y=ranginte*0.001d0
        call bpolyeder1(kmag,x,y,z,bx,by,bz)
        bxi(4)=bx/3.0d0
        byi(4)=by/3.0d0
        bzi(4)=bz/3.0d0

        y=ranginti*0.001d0-ds
        do i=2,nstepint-1,2
          y=y+ds2
          call bpolyeder1(kmag,x,y,z,bx,by,bz)
          bxi(2)=bxi(2)+bx
          byi(2)=byi(2)+by
          bzi(2)=bzi(2)+bz
        enddo

        y=ranginti*0.001d0
        do i=3,nstepint-2,2
          y=y+ds2
          call bpolyeder1(kmag,x,y,z,bx,by,bz)
          bxi(3)=bxi(3)+bx
          byi(3)=byi(3)+by
          bzi(3)=bzi(3)+bz
        enddo

        bxint=
     &    (bxi(1)/3.0d0+4.0d0/3.0d0*bxi(2)+2.0d0/3.0d0*bxi(3)+bxi(4)/3.0d0)*ds
        byint=
     &    (byi(1)/3.0d0+4.0d0/3.0d0*byi(2)+2.0d0/3.0d0*byi(3)+bxi(4)/3.0d0)*ds
        bzint=
     &    (bzi(1)/3.0d0+4.0d0/3.0d0*bzi(2)+2.0d0/3.0d0*bzi(3)+bxi(4)/3.0d0)*ds

      else if (vny.lt.tiny2.and.vnz.lt.tiny2) then

        y=yint
        z=zint

        x=ranginti*0.001d0
        call bpolyeder1(kmag,x,y,z,bx,by,bz)
        bxi(1)=bx/3.0d0
        byi(1)=by/3.0d0
        bzi(1)=bz/3.0d0

        x=ranginte*0.001d0
        call bpolyeder1(kmag,x,y,z,bx,by,bz)
        bxi(4)=bx/3.0d0
        byi(4)=by/3.0d0
        bzi(4)=bz/3.0d0

        x=ranginti*0.001d0-ds
        do i=2,nstepint-1,2
          x=x+ds2
          call bpolyeder1(kmag,x,y,z,bx,by,bz)
          bxi(2)=bxi(2)+bx
          byi(2)=byi(2)+by
          bzi(2)=bzi(2)+bz
        enddo

        x=ranginti*0.001d0
        do i=3,nstepint-2,2
          x=x+ds2
          call bpolyeder1(kmag,x,y,z,bx,by,bz)
          bxi(3)=bxi(3)+bx
          byi(3)=byi(3)+by
          bzi(3)=bzi(3)+bz
        enddo

        bxint=
     &    (bxi(1)/3.0d0+4.0d0/3.0d0*bxi(2)+2.0d0/3.0d0*bxi(3)+bxi(4)/3.0d0)*ds
        byint=
     &    (byi(1)/3.0d0+4.0d0/3.0d0*byi(2)+2.0d0/3.0d0*byi(3)+bxi(4)/3.0d0)*ds
        bzint=
     &    (bzi(1)/3.0d0+4.0d0/3.0d0*bzi(2)+2.0d0/3.0d0*bzi(3)+bxi(4)/3.0d0)*ds

      else ! vnxrot,vnyrot,vnzot parallel to an axis

          print*,'*** Warning in BPOLYINTNUM:'
          print*,
     & '*** line of integration not parallel to axis of coord.-system!'
          print*,'*** not yet implemented!'
          print*,'*** zero result returned'
          return

        endif ! vnxrot,vnyrot,vnzrot parallel to an axis

      return
      end
