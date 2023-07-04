*CMZ :  4.00/15 19/03/2022  09.37.41  by  Michael Scheer
*CMZ :  4.00/13 07/11/2021  17.42.14  by  Michael Scheer
*CMZ :  4.00/07 09/07/2020  12.43.50  by  Michael Scheer
*CMZ :  3.03/02 18/03/2016  16.15.39  by  Michael Scheer
*CMZ :  3.02/05 13/04/2015  11.55.17  by  Michael Scheer
*CMZ :  3.01/02 31/07/2013  14.44.36  by  Michael Scheer
*CMZ :  3.01/01 31/07/2013  12.20.54  by  Michael Scheer
*CMZ :  3.01/00 18/07/2013  13.37.20  by  Michael Scheer
*-- Author :    Michael Scheer   12/07/2013
      subroutine bue(xin,yin,zin,bxout,byout,bzout,axout,ayout,azout,imag)
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


*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,primkin.
      include 'primkin.cmn'
*KEEP,mgsqc.
      include 'mgsqc.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,halbasy.
      include 'halbasy.cmn'
*KEEP,ustep.
      include 'ustep.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      double precision xin,yin,zin,bxout,byout,bzout,axout,ayout,azout,
     &  halfshift,bx,by,bz,ax,ay,az,park,beff,consttap,dbytap,wlen1,dener,pdum,
     &  b0old,pkold,xlold,ylold,zlold,fasymold,ahwpold,
     &  hhold,xcenold,b0h,b0v,rhv,b0eff,
     &  phi,sphi(nmgsqp),cphi(nmgsqp),xphi,zphi,xx,zz

      integer imag,iwarn,ical,i,iahwfouro,nhold

      save ical
      data iwarn/0/
      data ical/0/

      bxout=0.0d0
      byout=0.0d0
      bzout=0.0d0

      axout=0.0d0
      ayout=0.0d0
      azout=0.0d0

      b0old=b0halbasy
      pkold=pkhalbasy
      xlold=xlhalbasy
      ylold=ylhalbasy
      zlold=zlhalbasy
      fasymold=fasym
      ahwpold=ahwpol
      iahwfouro=iahwfour
      nhold=nhhalbasy
      hhold=hhalbasy
      xcenold=xcenhal

      fasym=2.0d0
      pkhalbasy=0.0d0
      park=pmag(1,imag)
      b0v=pmag(2,imag)
      b0h=pmag(3,imag)
      zlhalbasy=pmag(6,imag)
      xlhalbasy=pmag(8,imag)
      nhhalbasy=pmag(9,imag)
      hhalbasy=pmag(10,imag)

      if (nhhalbasy.ne.0.and.hhalbasy.ne.0.0d0) then
        if (hhalbasy.eq.-9999.0d0) then
          if (ifreq2p.eq.1) then
            hhalbasy=freqlow
          else
            hhalbasy=(freqlow+freqhig)/2.0d0
          endif
        else if (hhalbasy.lt.0.0d0) then
          hhalbasy=-wtoe1/hhalbasy
        endif
        WLEN1=wtoe1/abs(hhalbasy/nhhalbasy)
        park=2.0d0*(wlen1/(zlhalbasy*1.0D9/2.0d0/dmygammap**2)-1.0d0)
        if (park.lt.0.0d0) then
          write(6,*)
     &      '*** Error in BUE:'
          write(6,*)
     &      'Inconsistent values for undulator of file magseq.in!'
          write(6,*)' '
          write(lungfo,*)
     &      '*** Error in BUE:'
          write(lungfo,*)
     &      'Inconsistent values for undulator of file magseq.in!'
          write(lungfo,*)' '
          stop
        endif
        park=sqrt(park)
      endif

      IF (park.NE.0.0) THEN

        B0EFF=park/(echarge1*zlhalbasy/(2.*PI1*EMASSKG1*CLIGHT1))

        if (b0h.eq.0.0d0.and.b0v.eq.0d0) then
          b0v=b0eff/sqrt(2.0d0)
          b0h=b0v
        else if (b0h.eq.0.0d0.and.b0v.ne.0d0) then
          b0v=b0v/abs(b0v)*b0eff
        else if (b0v.eq.0.0d0.and.b0h.ne.0d0) then
          b0h=b0h/abs(b0h)*b0eff
        else

          rhv=b0h/b0v

          b0h=b0eff/sqrt(1.0d0+1.0d0/rhv**2)*b0h/abs(b0h)
          b0v=b0h/rhv

        endif

        ! To avoid repetitions of initializations
        pmag(1,imag)=0.0d0
        pmag(2,imag)=b0v
        pmag(3,imag)=b0h
        pmag(9:10,imag)=0.0d0

      ENDIF

      if (b0v.eq.0.0d0.and.b0h.eq.0.0d0) return

      zkhalbasy=twopi1/zlhalbasy

      if(xlhalbasy.ne.0.0d0) then
        xkhalbasy=twopi1/xlhalbasy
        YKHALBASY=DSQRT(ZKHALBASY**2+XKHALBASY**2)
        ylhalbasy=twopi1/ykhalbasy
      else
        xkhalbasy=0.0d0
        YKHALBASY=ZKHALBASY
        YLHALBASY=ZLHALBASY
      endif

      ahwpol=(pmag(7,imag)-1.0d0)*2.0d0+1.0d0
      halfshift=pmag(4,imag)/2.0d0*zlhalbasy
      xcenhal=pmag(5,imag)-halfshift

c{ taper to adapt field strength to energy loss
      dbytap=0.0d0

      if (pmag(11,imag).ne.0.0d0) then

        if (ieneloss.eq.0.and.iwarn.eq.0) then
          write(lungfo,*)' '
          write(lungfo,*)'*** Warning in BUE: Taper found, but IENELOSS is zero!'
          write(lungfo,*)'*** Note: This factor scales the whole device. It is ment for a series of undulators ***'
          write(lungfo,*)' '
          write(6,*)' '
          write(6,*)'*** Warning in BUE: Taper found, but IENELOSS is zero!'
          write(6,*)'*** Note: This factor scales the whole device. It is ment for a series of undulators ***'
          write(6,*)' '
          iwarn=1
        endif

        beff=sqrt(b0h**2+b0v**2)
        park=beff*echarge1*zlhalbasy/(2.0d0*pi1*emasskg1*clight1)
        wlen1=(1+park**2/2.0d0)/2.0d0/dmygammaP**2*zlhalbasy*1.0d9
        consttap=dmyenergyP**2/(1.0d0+park**2/2.0d0)
        pdum=0.5d-9/pi1*cgam1*dmycur*(clight1*emassg1)**2
        dener=-pdum*beff**2/2.0d0
     &    *(ahwpol*zlhalbasy/2.0d0)
     &    /dmycur/1.0d9
     &    *dmygammaP**2 !loss per undu in GeV
        if (dmygammaP.ne.0.0d0) then
          dener=dmyenergyP*(gammaustep-dmygammaP)/dmygammaP !actual loss
        endif
        if (park.ne.0.0d0) then
          dbytap=pmag(11,imag)*
     &      2.0d0*dmyenergyP*dener/park**2/consttap !correction factor for taper
        else
          dbytap=0.0d0
        endif
      endif
      b0halbasy=b0v*(1.0d0+dbytap)
c} taper

      if (ical.eq.0) then
        do i=1,mmag
          phi=pmag(13,i)
          sphi(i)=sin(phi)
          cphi(i)=cos(phi)
        enddo
        ical=1
      endif

      if (b0v.ne.0.0d0) then
        pkhalbasy=0.0d0
        phi=pmag(13,imag)
        xx=xin-xcenhal
        zz=zin-pmag(12,imag)
        if (phi.ne.0.0d0) then
          xphi= cphi(imag)*xx+sphi(imag)*zz
          zphi=-sphi(imag)*xx+cphi(imag)*zz
        else
          xphi=xx
          zphi=zz
        endif
        call bhalbasy2(xphi,yin,zphi,bx,by,bz,ax,ay,az)
        bxout=bxout+bx
        byout=byout+by
        bzout=bzout+bz
        axout=axout+ax
        ayout=ayout+ay
        azout=azout+az
      endif
      if (b0h.ne.0.0d0) then
        b0halbasy=b0h*(1.0d0+dbytap)
        xcenhal=pmag(5,imag)+halfshift
        phi=pmag(13,imag)
        xx=xin-xcenhal
        zz=zin-pmag(12,imag)
        if (phi.ne.0.0d0) then
          xphi= cphi(imag)*xx+sphi(imag)*zz
          zphi=-sphi(imag)*xx+cphi(imag)*zz
        else
          xphi=xx
          zphi=zz
        endif
        call bhalbasy2(xphi,yin,zphi,bx,by,bz,ax,ay,az)
        bxout=bxout-bx
        byout=byout+bz
        bzout=bzout-by
        axout=axout+ax
        ayout=ayout+ay
        azout=azout+az
      endif

      b0halbasy=b0old
      pkhalbasy=pkold
      xlhalbasy=xlold
      ylhalbasy=ylold
      zlhalbasy=zlold
      fasym=fasymold
      ahwpol=ahwpold
      iahwfour=iahwfouro
      nhhalbasy=nhold
      hhalbasy=hhold
      xcenhal=xcenold

      return
      end
