*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.66/02 26/10/2009  14.33.39  by  Michael Scheer
*CMZ :  2.65/01 08/10/2009  09.58.11  by  Michael Scheer
*CMZ :  2.65/00 18/09/2009  08.30.19  by  Michael Scheer
*CMZ :  2.64/06 15/09/2009  09.52.14  by  Michael Scheer
*CMZ :  2.64/05 14/09/2009  10.06.11  by  Michael Scheer
*-- Author :    Michael Scheer   01/09/2009
      subroutine cyltocart(isour)
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

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,spectf90.
      include 'spectf90.cmn'
*KEND.

c
c Alte Version 154 wieder  hervorgeholt, da sie besser mit MPINR=0 uebereinstimmt.
c
c

c fill AFREQ and SPECPOW by interpolating AFREQRPHI and SPECPOWRPHI

      double precision, dimension(:,:,:), allocatable :: fr,coefr
      double precision, dimension(:,:), allocatable :: phas
      double precision, dimension(:), allocatable :: w1,w2,w3,w4,
     &  fp,coefp,obsvedge,frw

      double complex af2,af3
c     &  ,expom

      double precision y,z,r,phi,pihalf,phiedge,yp12,
     &  phiy,phiz,phiy1,phiz1,
     &  dist0,dist02,ddist,h2,cenxexi,dvlen,
     &  dphase,dphi,pi2phaser2,pi2phaser3,pi2phasephi2,pi2phasephi3,
     &  af2r,af2i,af3r,af3i,
     &  af2ro,af2io,af3ro,af3io

      integer ixy,iphi,ir,ifreq,iphas,iobrp,istat,medge,
     &  ifirst,ilast,iedge,isour,jphi,nphi,ical

      data ical/0/

      if (ical.eq.0) then

        medge=3
        yp12=9999.0d0

        allocate(phas(4,nobsv))

        allocate(frw(nobsvr))
        allocate(fr(nobsvr,4,nobsvphi))
        allocate(coefr(nobsvr,4,nobsvphi))

        allocate(fp(nobsvphi+2*medge+1))
        allocate(obsvedge(nobsvphi+2*medge+1))
        allocate(coefp(nobsvphi+2*medge+1))

        allocate(w1(max(nobsvr,nobsvphi+2*medge+1)))
        allocate(w2(max(nobsvr,nobsvphi+2*medge+1)))
        allocate(w3(max(nobsvr,nobsvphi+2*medge+1)))
        allocate(w4(max(nobsvr,nobsvphi+2*medge+1)))

        pihalf=pi1/2.0d0

        ifirst=medge+1

        if (iquadphi.eq.0) then
          ilast=ifirst+nobsvphi
          obsvedge(ilast)=obsvphi(nobsvphi)+obsvdphi
        else
          ilast=ifirst+nobsvphi-1
        endif

        do iphi=1,nobsvphi
          obsvedge(ifirst+iphi-1)=obsvphi(iphi)
        enddo

        do iedge=1,medge
          obsvedge(ilast+iedge)=obsvedge(ilast+iedge-1)+obsvdphi
        enddo

        do iedge=1,medge
          obsvedge(ifirst-iedge)=obsvedge(ifirst-iedge+1)-obsvdphi
        enddo

        ical=1

      endif

        cenxexi=(min(sourceeo(1,1,isour),xiend)
     &    +max(sourceao(1,1,isour),xianf))/2.d0
        dvlen=
     &    min(sourceeo(1,1,isour),xiend)-
     &    max(sourceao(1,1,isour),xianf)
        if (iampli.lt.0) then
          dvlen=dvlen*(-iampli)
        endif
        dist0=pincen(1)-cenxexi
        dist02=dist0**2

      do ifreq=1,nfreq

        iobrp=0
        ixy=1

        pi2phasephi2=0.0d0
        pi2phasephi3=0.0d0
        do iphi=1,nobsvphi

          phiy1=-twopi1
          phiz1=-twopi1

          pi2phaser2=0.0d0
          pi2phaser3=0.0d0
          do ir=1,nobsvr

            iobrp=iobrp+1
            ifrob=ifreq+nfreq*(iobrp-1)

            h2=(obsvr(ir)/dist0)**2
            if (h2.lt.0.01) then
              ddist=dist0*(h2/2.0d0-h2**2/8.0d0)
            else
              ddist=dist0*(sqrt(h2)-1.0d0)
            endif

            dphase=ddist/freq(ifreq)*wtoe1*1.0d9*twopi1

            af2=afreqrphi(2,IFROB)
            af3=afreqrphi(3,IFROB)
            af2r=dreal(af2)
            af2i=dimag(af2)
            af3r=dreal(af3)
            af3i=dimag(af3)

            phiy=atan2(af2i,af2r)
            phiz=atan2(af3i,af3r)

            fr(ir,1,iphi)=abs(af2)
            fr(ir,2,iphi)=phiy-dphase
            fr(ir,3,iphi)=abs(af3)
            fr(ir,4,iphi)=phiz-dphase

c            print*,ir,iphi,phiz+dphase,phiz-dphase

            phiy1=phiy
            phiz1=phiz
            af2ro=af2r
            af2io=af2i
            af3ro=af3r
            af3io=af3i

          enddo !nobsvr

          do iphas=2,4,2
            frw(1)=fr(1,iphas,iphi)
            do ir=2,nobsvr
              frw(ir)=fr(ir,iphas,iphi)
              phiy1=frw(ir-1)
              phiy=frw(ir)
              dphi=phiy-phiy1

              if (dphi.gt.pi1) then
                dphi=dphi-twopi1
              else if (dphi.lt.-pi1) then
                dphi=dphi+twopi1
              endif
              fr(ir,iphas,iphi)=fr(ir-1,iphas,iphi)+dphi

            enddo
          enddo

          do iphas=1,4
            call util_spline_coef(obsvr,
     &        fr(1,iphas,iphi),nobsvr,yp12,yp12,coefr(1,iphas,iphi),
     &        w1,w2,w3,w4)
          enddo

        enddo !nobsvphi

        do iphas=1,4

          do ixy=1,nobsv

            y=obsv(2,ixy)
            z=obsv(3,ixy)
            r=sqrt(y**2+z**2)

            if (iquadphi.eq.0) then
              phi=atan2(y,z)
            else
              if (z.ne.0.0d0) then
                phi=atan(y/z)
              else
                phi=pi1/2.0d0
              endif
              phi=abs(phi)
            endif

            if (phi.lt.0.0d0) phi=phi+twopi1

            do iphi=1,nobsvphi

              call util_spline_inter_status(
     &          obsvr,fr(1,iphas,iphi),coefr(1,iphas,iphi),
     &          nobsvr,r,fp(iphi+medge),0,
     &          istat)

              if (istat.ne.0) then
                stop
     &            '*** Error: Bad return from util_spline_inter_status in CYLTOCART'
              endif

            enddo !iphi

            if (iquadphi.eq.0) then

              nphi=nobsvphi+1+2*medge
              ifirst=medge+1
              ilast=ifirst+nobsvphi

              fp(ilast)=fp(ifirst) !Periode vervollständigen

              do iedge=1,medge
                fp(ilast+iedge)=fp(ifirst+iedge)
              enddo

              do iedge=1,medge
                fp(ifirst-iedge)=fp(ilast-iedge)
              enddo

            else !if (iquadphi.eq.0)

              nphi=nobsvphi+2*medge
              ifirst=medge+1
              ilast=ifirst+nobsvphi-1

              do iedge=1,medge
                phiedge=mod(pihalf+iedge*obsvdphi,twopi1)
                if (phiedge.gt.3.0d0*pihalf) then
                  phiedge=twopi1-phiedge
                else if (phiedge.gt.pi1) then
                  phiedge=phiedge-pi1
                else if (phiedge.gt.pihalf) then
                  phiedge=pi1-phiedge
                endif
                jphi=ifirst+nint(phiedge/obsvdphi)
                fp(ilast+iedge)=fp(jphi)
              enddo

              do iedge=1,medge
                phiedge=mod(iedge*obsvdphi,twopi1)
                if (phiedge.gt.3.0d0*pihalf) then
                  phiedge=twopi1-phiedge
                else if (phiedge.gt.pi1) then
                  phiedge=phiedge-pi1
                else if (phiedge.gt.pihalf) then
                  phiedge=pi1-phiedge
                endif
                jphi=ifirst+nint(phiedge/obsvdphi)
                fp(ifirst-iedge)=fp(jphi)
              enddo

            endif !iquadphi

            call util_spline_coef(obsvedge,fp,nphi,
     &        yp12,yp12,coefp,
     &        w1,w2,w3,w4)

            call util_spline_inter_status(
     &        obsvedge,fp,coefp,nphi,phi,phas(iphas,ixy),0,istat)

            if (iphas.eq.4) then

              h2=(r/dist0)**2

              if (h2.lt.0.01) then
                ddist=dist0*(h2/2.0d0-h2**2/8.0d0)
              else
                ddist=dist0*(sqrt(h2)-1.0d0)
              endif

              dphase=ddist/freq(ifreq)*wtoe1*1.0d9*twopi1

              afreq(2,ifreq+nfreq*(ixy-1))=
     &          phas(1,ixy)*
     &          cmplx(cos(phas(2,ixy)+dphase),sin(phas(2,ixy)+dphase))

              afreq(3,ifreq+nfreq*(ixy-1))=
     &          phas(3,ixy)*
     &          cmplx(cos(phas(4,ixy)+dphase),sin(phas(4,ixy)+dphase))

            endif

            if (istat.ne.0) then
              stop
     &          '*** Error: Bad return from util_spline_inter_status in CYLTOCART'
            endif

          enddo !ixy

        enddo !iphas

      enddo !ifreq

c power

      iphas=1
      do ixy=1,nobsv

        y=obsv(2,ixy)
        z=obsv(3,ixy)
        r=sqrt(y**2+z**2)

        if (iquadphi.eq.0) then
          phi=atan2(y,z)
        else
          if (z.ne.0.0d0) then
            phi=atan(y/z)
          else
            phi=pi1/2.0d0
          endif
          phi=abs(phi)
        endif

        if (phi.lt.0.0d0) phi=phi+twopi1

        iobrp=0
        do iphi=1,nobsvphi

          if (ixy.eq.1) then

            do ir=1,nobsvr
              iobrp=iobrp+1
              ifrob=isour+nsource*(iobrp-1)
              fr(ir,iphas,iphi)=specpowrphi(ifrob)
            enddo

            call util_spline_coef(obsvr,fr(1,iphas,iphi),nobsvr,yp12,yp12,
     &        coefr(1,iphas,iphi),
     &        w1,w2,w3,w4)

          endif !ixy.eq.1

          call util_spline_inter_status(
     &      obsvr,fr(1,iphas,iphi),coefr(1,iphas,iphi),nobsvr,r,
     &      fp(iphi+medge),0,istat)

          if (istat.ne.0) then
            stop
     &        '*** Error: Bad return from util_spline_inter_status in CYLTOCART'
          endif

        enddo !iphi

        if (iquadphi.eq.0) then

          ifirst=medge+1
          ilast=ifirst+nobsvphi
          fp(ilast)=fp(ifirst) !Periode vervollständigen

          do iedge=1,medge
            fp(ilast+iedge)=fp(ifirst+iedge)
          enddo

          do iedge=1,medge
            fp(ifirst-iedge)=fp(ilast-iedge)
          enddo

        else !if (iquadphi.eq.0)

          ifirst=medge+1
          ilast=ifirst+nobsvphi-1

          do iedge=1,medge
            phiedge=mod(pihalf+iedge*obsvdphi,twopi1)
            if (phiedge.gt.3.0d0*pihalf) then
              phiedge=twopi1-phiedge
            else if (phiedge.gt.pi1) then
              phiedge=phiedge-pi1
            else if (phiedge.gt.pihalf) then
              phiedge=pi1-phiedge
            endif
            jphi=ifirst+nint(phiedge/obsvdphi)
            fp(ilast+iedge)=fp(jphi)
          enddo

          do iedge=1,medge
            phiedge=mod(iedge*obsvdphi,twopi1)
            if (phiedge.gt.3.0d0*pihalf) then
              phiedge=twopi1-phiedge
            else if (phiedge.gt.pi1) then
              phiedge=phiedge-pi1
            else if (phiedge.gt.pihalf) then
              phiedge=pi1-phiedge
            endif
            jphi=ifirst+nint(phiedge/obsvdphi)
            fp(ifirst-iedge)=fp(jphi)
          enddo

        endif !iquadphi

        call util_spline_coef(obsvedge,fp,nphi,
     &    yp12,yp12,coefp,
     &    w1,w2,w3,w4)

        call util_spline_inter_status(
     &    obsvedge,fp,coefp,nphi,phi,
     &    specpow(isour+nsource*(ixy-1)),0,istat)

        if (istat.ne.0) then
          stop
     &      '*** Error: Bad return from util_spline_inter_status in CYLTOCART'
        endif

      enddo !ixy

      if (isour.eq.nsource) then
        deallocate(phas)
        deallocate(fr)
        deallocate(fp)
        deallocate(obsvedge)
        deallocate(coefr)
        deallocate(coefp)

        deallocate(w1)
        deallocate(w2)
        deallocate(w3)
        deallocate(w4)
      endif

      return
      end
