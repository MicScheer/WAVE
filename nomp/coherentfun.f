*CMZ :  2.66/12 25/10/2012  15.10.37  by  Michael Scheer
*-- Author :    Michael Scheer   06/05/2010
      function coherentfun(x)
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

      use bunchmod

      double precision phi,fun,sig,wlen,emod,dpsidp,p0,r56,bessel,freqr,wscale,
     &  dtphase
      real coherentfun,x
      integer ical,icoef,nsigs,iwarn,ifail,n

      data ical/1/
      data nsigs/10/

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cohfun.
      include 'cohfun.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,trackf90.
      include 'trackf90.cmn'
*KEND.

      if (ical.eq.1) then

        if (ibunch.ne.0) then
          p0=bunchp0
          r56=bunchr56
        else if (iampli.ne.0) then
          p0=ampbunchp0
          r56=ampbunchr56
        else
            write(lungfo,*)
     &        '*** Error in COHERENTFUN: IBUNCH nor IAMPLI .ne. 0'
            write(6,*)
     &        '*** Error in COHERENTFUN: IBUNCH nor IAMPLI .ne. 0'
            stop '*** Program WAVE aborted ***'
        endif

        sig=cfxlenfou/2.0d0/nsigs

        if (p0.le.0.0d0) p0=1.0d-30
        if (r56.le.0.0d0) r56=1.0d-30

        dtphase=(wtra2i+(1.d0/(dmygamma*dmygamma))
     &    *(xstop-xstart)/2.d0)/clight1
        freqr=2.d0*pi1/dtphase*hbarev1

        if (ampfreq.eq.-9999.d0) then
          ampfreq=freqr
        endif  !(imampli.gt.0.and.ampfreq.eq.-9999.d0)

        wlen=wtoe1/ampfreq*1.0d-9

        ical=1
      endif

      if (modecf.eq.1) then

        coherentfun=exp(-(x/sig)**2/2.0d0)/sqrttwopi1

      else if (modecf.eq.2) then

          dpsidp=twopi1/wlen*r56/dmyenergy
          wscale=0.0d0

          do n=1,ncoefcf
            emod=exp(-(n*espread*dmyenergy*dpsidp)**2/2.0d0)
            call util_bessel(n,n*p0*dpsidp,bessel,ifail)
            if (ifail.ne.0.and.iwarn.eq.0) then
              write(lungfo,*)'*** Warning: Error in COHERENTFUN from call to UTIL_BESSEL ***'
              write(lungfo,*)'*** Warning: further warnings will be suppressed ***'
              write(6,*)'*** Warning: Error in COHERENTFUN from call to UTIL_BESSEL ***'
              write(6,*)'*** Warning: further warnings will be suppressed ***'
              iwarn=1
            endif
            wscale=wscale+2.0d0*emod*bessel
          enddo

          if (wscale.gt.0.0d0) then
            wscale=1.0d0/wscale
          else
            wscale=1.0d0
          endif

          coherentfun=wscale

      else if (modecf.eq.4) then

        fun=ccf(0)
        phi=twopi1*x/cfxlenfou

        do icoef=1,ncoefcf
          fun=fun+ccf(icoef)*cos(icoef*phi)
        enddo
        coherentfun=fun

      else
        stop
     & '*** Error in COHERENTFUN: Undefined mode, please check IUBUNCH or IMAMPLI'
      endif

      return
      end
