*CMZ :  3.05/00 25/04/2018  13.09.51  by  Michael Scheer
*CMZ :  2.52/06 19/08/2004  15.34.24  by  Michael Scheer
*CMZ :  2.52/05 16/08/2004  13.37.21  by  Michael Scheer
*CMZ :  1.02/01 09/08/2004  14.45.40  by  Michael Scheer
*CMZ :  1.02/00 28/07/2004  16.57.47  by  Michael Scheer
*-- Author :    Michael Scheer   27/07/2004
      subroutine bpolyint(kmag,xint,yint,zint,
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

c literature: Elleaume, Chubar, Chavanne PAC97
c             Computing 3D Magnetic Fields form Insertion Devices

c calculates first magnetic integral for line defined by point
c      (xint,yint,zint) and vector (vxint,vyint,vzint)
c restrictions:
c     rectangular magnets only
c     line of integrations parallel to axis of block and block aligned to
c     axis of coordinate system (paper gives formulas for general
c     case of integration direction (not coded here))

c xint,yint,zint given in meter
c bxint,byint,bzint given in Tm

*KEEP,bwpolyederf90u.
      include 'bwpolyederf90u.cmn'
*KEND.

      implicit none

      double precision xint,yint,zint,vxint,vyint,vzint,vnx,vny,vnz,vn,
     &  bxint,byint,bzint,
     &  vnxrot,vnyrot,vnzrot,
     &  vnxrota,vnyrota,vnzrota,
     &  xr(2),yr(2),zr(2),ts(3,3),
     &  tsinv(3,3),w(3),x12,x22,z12,z22,y12,y22

      double precision g(3,3),vmagrot(3),vmaglab(3),bint1(3),xxrot,yyrot,zzrot,
     &  xxint,yyint,zzint,tiny2,pi2inv,pi4inv

      parameter (pi2inv=0.159154943091895d0)
      parameter (pi4inv=0.0795774715459477d0)

      integer ical
      integer kmag,i,k

      data ical/0/

      if (ical.eq.0) then
        tiny2=tiny*tiny
        ical=1
      endif !ical

      xxint=xint*1000.0d0
      yyint=yint*1000.0d0
      zzint=zint*1000.0d0

      vn=sqrt(vxint*vxint+vyint*vyint+vzint*vzint)

      vnx=vxint/vn
      vny=vyint/vn
      vnz=vzint/vn

      bxint=0.d0
      byint=0.d0
      bzint=0.d0

      if (bpebc(7,kmag).eq.0.0d0) return

      if(bpebc(8,kmag).ne.-6) then !not rectangular magnet

        print*,'*** Warning in BPOLYINT: Non-rectangular magnet!'
        print*,'***                      not yet implemented!'
        print*,'***                      zero result returned'
        return

      else !bpebc(8,kmag).ne.-6

        vmaglab(1)=bpebc(4,kmag)
        vmaglab(2)=bpebc(5,kmag)
        vmaglab(3)=bpebc(6,kmag)

c transform everything to the nz=(0,0,1) system and rotate it parallel to x-axis

        do i=1,3
          do k=1,3
            ts(i,k)=bpetm(i,k,1,kmag)
            tsinv(i,k)=bpetm(i,k+3,1,kmag)
          enddo
        enddo

        xxrot=ts(1,1)*xxint+ts(1,2)*yyint+ts(1,3)*zzint
        yyrot=ts(2,1)*xxint+ts(2,2)*yyint+ts(2,3)*zzint
        zzrot=ts(3,1)*xxint+ts(3,2)*yyint+ts(3,3)*zzint

        vnxrot=ts(1,1)*vnx+ts(1,2)*vny+ts(1,3)*vnz
        vnyrot=ts(2,1)*vnx+ts(2,2)*vny+ts(2,3)*vnz
        vnzrot=ts(3,1)*vnx+ts(3,2)*vny+ts(3,3)*vnz

        vnxrota=abs(vnxrot)
        vnyrota=abs(vnyrot)
        vnzrota=abs(vnzrot)

        vmagrot(1)=
     &    ts(1,1)*vmaglab(1)+ts(1,2)*vmaglab(2)+ts(1,3)*vmaglab(3)
        vmagrot(2)=
     &    ts(2,1)*vmaglab(1)+ts(2,2)*vmaglab(2)+ts(2,3)*vmaglab(3)
        vmagrot(3)=
     &    ts(3,1)*vmaglab(1)+ts(3,2)*vmaglab(2)+ts(3,3)*vmaglab(3)

c dimensions of magnet

        w(1)=bperot(1,1,1,kmag)-bperot(1,2,1,kmag)
        w(2)=bperot(2,1,1,kmag)-bperot(2,3,1,kmag)
        w(3)=bperot(3,1,1,kmag)-bperot(3,1,3,kmag)

c distances from considered point to corners of magnet

        xr(1)=bperot(1,1,1,kmag)-xxrot
        xr(2)=bperot(1,2,1,kmag)-xxrot
        yr(1)=bperot(2,1,1,kmag)-yyrot
        yr(2)=bperot(2,3,1,kmag)-yyrot

        zr(1)=bperot(3,1,1,kmag)-zzrot
        zr(2)=bperot(3,1,3,kmag)-zzrot

        g=0.0d0

        if (vnxrota.lt.tiny2.and.vnyrota.lt.tiny2) then

          if (abs(xr(1)).lt.tiny2) then
            print*,'Warning in BPOLYINT: xr(1) too small'
            print*,'magnet:',kmag
            print*,'xr(1)',xr(1)
            if (xr(1).lt.0.0d0) then
              print*,'Set to',-tiny2
              xr(1)=-tiny2
            else
              print*,'Set to',tiny2
              xr(1)=tiny2
            endif
          endif


          if (abs(xr(2)).lt.tiny2) then
            print*,'Warning in BPOLYINT: xr(2) too small'
            print*,'magnet:',kmag
            print*,'xr(2)',xr(2)
            if (xr(2).lt.0.0d0) then
              print*,'Set to',-tiny2
              xr(2)=-tiny2
            else
              print*,'Set to',tiny2
              xr(2)=tiny2
            endif
          endif

          if (abs(yr(1)).lt.tiny2) then
            print*,'Warning in BPOLYINT: yr(1) too small'
            print*,'magnet:',kmag
            print*,'yr(1)',yr(1)
            if (yr(1).lt.0.0d0) then
              print*,'Set to',-tiny2
              yr(1)=-tiny2
            else
              print*,'Set to',tiny2
              yr(1)=tiny2
            endif
          endif


          if (abs(yr(2)).lt.tiny2) then
            print*,'Warning in BPOLYINT: yr(2) too small'
            print*,'magnet:',kmag
            print*,'yr(2)',yr(2)
            if (yr(2).lt.0.0d0) then
              print*,'Set to',-tiny2
              yr(2)=-tiny2
            else
              print*,'Set to',tiny2
              yr(2)=tiny2
            endif
          endif

          do i=1,2
            do k=1,2
              g(1,1)=g(1,1)+(-1)**(i+k)*atan(xr(i)/yr(k))
              g(2,2)=g(2,2)+(-1)**(i+k)*atan(yr(k)/xr(i))
            enddo !k
          enddo !i

          x12=xr(1)*xr(1)
          x22=xr(2)*xr(2)
          y12=yr(1)*yr(1)
          y22=yr(2)*yr(2)

          g(1,2)=w(3)*pi4inv*log((x12+y22)*(x22+y12)/((x12+y12)*(x22+y22)))
     &          *0.001d0 !Tmm -> Tm
          g(2,1)=g(1,2)

          g(1,1)=g(1,1)*w(3)*pi2inv
     &          *0.001d0 !Tmm -> Tm
          g(2,2)=g(2,2)*w(3)*pi2inv
     &          *0.001d0 !Tmm -> Tm

        else if (vnxrota.lt.tiny2.and.vnzrota.lt.tiny2) then

          if (abs(xr(1)).lt.tiny2) then
            print*,'Warning in BPOLYINT: xr(1) too small'
            print*,'magnet:',kmag
            print*,'xr(1)',xr(1)
            if (xr(1).lt.0.0d0) then
              print*,'Set to',-tiny2
              xr(1)=-tiny2
            else
              print*,'Set to',tiny2
              xr(1)=tiny2
            endif
          endif

          if (abs(xr(2)).lt.tiny2) then
            print*,'Warning in BPOLYINT: xr(2) too small'
            print*,'magnet:',kmag
            print*,'xr(2)',xr(2)
            if (xr(2).lt.0.0d0) then
              print*,'Set to',-tiny2
              xr(2)=-tiny2
            else
              print*,'Set to',tiny2
              xr(2)=tiny2
            endif
          endif

          if (abs(zr(1)).lt.tiny2) then
            print*,'Warning in BPOLYINT: zr(1) too small'
            print*,'magnet:',kmag
            print*,'zr(1)',zr(1)
            if (zr(1).lt.0.0d0) then
              print*,'Set to',-tiny2
              zr(1)=-tiny2
            else
              print*,'Set to',tiny2
              zr(1)=tiny2
            endif
          endif

          if (abs(zr(2)).lt.tiny2) then
            print*,'Warning in BPOLYINT: zr(2) too small'
            print*,'magnet:',kmag
            print*,'zr(2)',zr(2)
            if (zr(2).lt.0.0d0) then
              print*,'Set to',-tiny2
              zr(2)=-tiny2
            else
              print*,'Set to',tiny2
              zr(2)=tiny2
            endif
          endif

          do i=1,2
            do k=1,2
              g(1,1)=g(1,1)+(-1)**(i+k)*atan(xr(i)/zr(k))
              g(3,3)=g(3,3)+(-1)**(i+k)*atan(zr(k)/xr(i))
            enddo !k
          enddo !i

          x12=xr(1)*xr(1)
          x22=xr(2)*xr(2)
          z12=zr(1)*zr(1)
          z22=zr(2)*zr(2)

          g(1,3)=w(2)*pi4inv*log((x12+z22)*(x22+z12)/((x12+z12)*(x22+z22)))
     &          *0.001d0 !Tmm -> Tm
          g(3,1)=g(1,3)

          g(1,1)=g(1,1)*w(2)*pi2inv
     &          *0.001d0 !Tmm -> Tm
          g(3,3)=g(3,3)*w(2)*pi2inv
     &          *0.001d0 !Tmm -> Tm

        else if (vnyrota.lt.tiny2.and.vnzrota.lt.tiny2) then

          if (abs(yr(1)).lt.tiny2) then
            print*,'Warning in BPOLYINT: yr(1) too small'
            print*,'magnet:',kmag
            print*,'yr(1)',yr(1)
            if (yr(1).lt.0.0d0) then
              print*,'Set to',-tiny2
              yr(1)=-tiny2
            else
              print*,'Set to',tiny2
              yr(1)=tiny2
            endif
          endif

          if (abs(yr(2)).lt.tiny2) then
            print*,'Warning in BPOLYINT: yr(2) too small'
            print*,'magnet:',kmag
            print*,'yr(2)',yr(2)
            if (yr(2).lt.0.0d0) then
              print*,'Set to',-tiny2
              yr(2)=-tiny2
            else
              print*,'Set to',tiny2
              yr(2)=tiny2
            endif
          endif

          if (abs(zr(1)).lt.tiny2) then
            print*,'Warning in BPOLYINT: zr(1) too small'
            print*,'magnet:',kmag
            print*,'zr(1)',zr(1)
            if (zr(1).lt.0.0d0) then
              print*,'Set to',-tiny2
              zr(1)=-tiny2
            else
              print*,'Set to',tiny2
              zr(1)=tiny2
            endif
          endif

          if (abs(zr(2)).lt.tiny2) then
            print*,'Warning in BPOLYINT: zr(2) too small'
            print*,'magnet:',kmag
            print*,'zr(2)',zr(2)
            if (zr(2).lt.0.0d0) then
              print*,'Set to',-tiny2
              zr(2)=-tiny2
            else
              print*,'Set to',tiny2
              zr(2)=tiny2
            endif
          endif

          do i=1,2
            do k=1,2
              g(2,2)=g(2,2)+(-1)**(i+k)*atan(yr(i)/zr(k))
              g(3,3)=g(3,3)+(-1)**(i+k)*atan(zr(k)/yr(i))
            enddo !k
          enddo !i

          y12=yr(1)*yr(1)
          y22=yr(2)*yr(2)
          z12=zr(1)*zr(1)
          z22=zr(2)*zr(2)

          g(2,3)=w(1)*pi4inv*log((y12+z22)*(y22+z12)/((y12+z12)*(y22+z22)))
     &          *0.001d0 !Tmm -> Tm
          g(3,2)=g(2,3)

          g(2,2)=g(2,2)*w(1)*pi2inv
     &          *0.001d0 !Tmm -> Tm
          g(3,3)=g(3,3)*w(1)*pi2inv
     &          *0.001d0 !Tmm -> Tm

        else ! vnxrot,vnyrot,vnzot parallel to an axis

          print*,'*** Warning in BPOLYINT:'
          print*,
     &      '*** line of integration not parallel to axis of coord.-system!'
          print*,'*** not yet implemented!'
          print*,'*** zero result returned'
          return

        endif ! vnxrot,vnyrot,vnzrot parallel to an axis

        bint1(1)=(g(1,1)*vmagrot(1)+g(1,2)*vmagrot(2)+g(1,3)*vmagrot(3))
        bint1(2)=(g(2,1)*vmagrot(1)+g(2,2)*vmagrot(2)+g(2,3)*vmagrot(3))
        bint1(3)=(g(3,1)*vmagrot(1)+g(3,2)*vmagrot(2)+g(3,3)*vmagrot(3))

        bxint=tsinv(1,1)*bint1(1)+tsinv(1,2)*bint1(2)+tsinv(1,3)*bint1(3)
        byint=tsinv(2,1)*bint1(1)+tsinv(2,2)*bint1(2)+tsinv(2,3)*bint1(3)
        bzint=tsinv(3,1)*bint1(1)+tsinv(3,2)*bint1(2)+tsinv(3,3)*bint1(3)

      endif !(bpebc(8,kmag).eq.1)

      return
      end
