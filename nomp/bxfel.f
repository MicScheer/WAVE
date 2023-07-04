*CMZ :  3.01/00 04/07/2013  08.32.02  by  Michael Scheer
*CMZ :  3.00/02 09/04/2013  15.30.14  by  Michael Scheer
*CMZ :  3.00/01 28/03/2013  10.01.42  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine bxfel(xin,y,z,bxout,byout,bzout,axout,ayout,azout)
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

c +PATCH,//WAVE/UFOR
c +DECK,bxfel.
      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,halbasy.
      include 'halbasy.cmn'
*KEND.

      double precision xin,x,y,z,bx,by,bz,ax,ay,az,
     &  bxout,byout,bzout,
     &  axout,ayout,azout,
     &  sectlen,consttap,dbytap,dener,ctaper,ctaper1,ctaper2,
     &  defl,zperlenm,zperlencm,
     &  bdefl,pdum,eharm1,park,wlen1,ahwpolundu

      integer ical,iundu

      data ical/0/

      if (ical.eq.0) then

        zperlenm=user(5)
        sectlen=user(4)
        zperlencm=zperlenm*100.0d0

        bdefl=ECHARGE1*zperlenm/(2.0d0*pi1*emasskg1*clight1)
        defl=bdefl*user(1)
        park=defl
        park=echarge1*dabs(user(1))*zperlenm/(2.*pi1*emasskg1*clight1)
        wlen1=(1+park**2/2.)/2./dmygamma**2*zlhalbasy*1.d9
        if (wlen1.ne.0.0d0) eharm1=wtoe1/wlen1

        consttap=dmyenergy**2/(1.0d0+defl**2/2.0d0)
        ahwpolundu=user(6)

        pdum=0.5d-9/pi1*cgam1*dmycur*(clight1*emassg1)**2
        dener=-pdum*user(1)**2/2.0d0
     &    *(ahwpolundu*zperlenm/2.0d0)
     &    /dmycur/1.0d9
     &    *dmygamma**2 !loss per undu in GeV

        if (defl.ne.0.0d0) then
          dbytap=user(3)*2.0d0*dmyenergy*dener/defl/consttap/bdefl !correction factor for taper
        else
          dbytap=0.0d0
        endif

        write(16,*)' '
        write(16,*)'     Subroutine bxfel called:'
        write(16,*)' '
        write(16,*)'     B0 [T] and K of undulator:             ',user(1), defl
        write(16,*)'     B0 of phaseshifter:                    ',user(2)
        write(16,*)'     Scaling of taper:                      ',user(3)
        write(16,*)'     Length of section [m]:                 ',user(4)
        write(16,*)'     Period-length of undulators [m]:       ',user(5)
        write(16,*)'     X-position of first undulator [m]:     ',user(7)
        write(16,*)'     X-position of second undulator [m]:    ',user(8)
        write(16,*)'     Period-length of phase-shifters [m]:   ',user(9)
        write(16,*)'     X-position of first phase-shifter [m]: ',user(10)
        write(16,*)'     X-position of second phase-shifter [m]:',user(11)
        write(16,*)' '

        if (user(3).ne.0.0d0.and.ieneloss.eq.0) then
          write(6,*)' '
          write(6,*)'*** Warning in BXFEL: Taper is not zero, but IENELOSS!'
          write(6,*)' '
          write(16,*)' '
          write(16,*)'*** Warning in BXFEL: Taper is not zero, but IENELOSS!'
          write(16,*)' '
        endif

        write(16,*)'     Energy loss per undulator [GeV]:',-dener

        if (user(1).ne.0.0) then
          write(16,*)'     Rel. taper for By:',dbytap/user(1)
        else
          write(16,*)'     Rel. taper for By: 0.0'
        endif

        ical=1
      endif !ical

      x=mod(xin,sectlen)

      bxout=0.0d0
      byout=0.0d0
      bzout=0.0d0
      axout=0.0d0
      ayout=0.0d0
      azout=0.0d0

      iundu=xin/(sectlen/2.0d0)
      if (user(1).ne.0.0d0) then
        ctaper=1.0d0+(dbytap*iundu)/user(1)
      else
        ctaper=1.0d0
      endif
      ctaper1=ctaper

c first undulator

      nhhalbasy=0
      pkhalbasy=0.0d0
      b0halbasy=user(1)*ctaper
      xlhalbasy=0.0d0
      xkhalbasy=0.0d0
      zlhalbasy=zperlenm
      zkhalbasy=twopi1/zlhalbasy
      ylhalbasy=zlhalbasy
      ykhalbasy=zkhalbasy
      ahwpol=ahwpolundu
      xcenhal=user(7)
      call bhalbasy(x,y,z,bx,by,bz,ax,ay,az)
      bxout=bxout+bx
      byout=byout+by
      bzout=bzout+bz
      axout=axout+ax
      ayout=ayout+ay
      azout=azout+az

c first phase-shifter
      nhhalbasy=0
      pkhalbasy=0.0d0
      b0halbasy=user(2)*ctaper
      xlhalbasy=0.0d0
      xkhalbasy=0.0d0
      zlhalbasy=user(9)
      zkhalbasy=twopi1/zlhalbasy
      ylhalbasy=zlhalbasy
      ykhalbasy=zkhalbasy
      ahwpol=1.0d0
      xcenhal=user(10)
      call bhalbasy(x,y,z,bx,by,bz,ax,ay,az)
      bxout=bxout+bx
      byout=byout+by
      bzout=bzout+bz
      axout=axout+ax
      ayout=ayout+ay
      azout=azout+az

      iundu=iundu+1
      if (user(1).ne.0.0d0) then
        ctaper=1.0d0+(dbytap*iundu)/user(1)
      else
        ctaper=1.0d0
      endif
      ctaper2=ctaper

c second undulator
      nhhalbasy=0
      pkhalbasy=0.0d0
      b0halbasy=user(1)*ctaper
      xlhalbasy=0.0d0
      xkhalbasy=0.0d0
      zlhalbasy=zperlenm
      zkhalbasy=twopi1/zlhalbasy
      ylhalbasy=zlhalbasy
      ykhalbasy=zkhalbasy
      ahwpol=ahwpolundu
      xcenhal=user(8)
      call bhalbasy(x,y,z,bx,by,bz,ax,ay,az)
      bxout=bxout+bx
      byout=byout+by
      bzout=bzout+bz
      axout=axout+ax
      ayout=ayout+ay
      azout=azout+az

c second phase-shifter
      nhhalbasy=0
      pkhalbasy=0.0d0
      b0halbasy=user(2)*ctaper
      xlhalbasy=0.0d0
      xkhalbasy=0.0d0
      zlhalbasy=user(9)
      zkhalbasy=twopi1/zlhalbasy
      ylhalbasy=zlhalbasy
      ykhalbasy=zkhalbasy
      ahwpol=1.0d0
      xcenhal=user(11)
      call bhalbasy(x,y,z,bx,by,bz,ax,ay,az)
      bxout=bxout+bx
      byout=byout+by
      bzout=bzout+bz
      axout=axout+ax
      ayout=ayout+ay
      azout=azout+az

c quadrupoles

      if (x.gt.sectlen/2.0d0) then
        ctaper=ctaper2
      else
        ctaper=ctaper1
      endif
      call bmagseq(x,y,z,bx,by,bz,ax,ay,az)
      bxout=bxout+bx*ctaper
      byout=byout+by*ctaper
      bzout=bzout+bz*ctaper
      axout=axout+ax*ctaper
      ayout=ayout+ay*ctaper
      azout=azout+az*ctaper

      return
      end
