*CMZ :  3.06/00 14/02/2019  21.01.16  by  Michael Scheer
*CMZ :  3.02/03 03/11/2014  14.55.22  by  Michael Scheer
*CMZ :  3.02/00 10/09/2014  09.38.31  by  Michael Scheer
*CMZ :  3.01/08 01/07/2014  11.04.35  by  Michael Scheer
*CMZ :  2.66/12 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.66/09 25/03/2010  13.07.39  by  Michael Scheer
*CMZ :  2.66/08 16/03/2010  13.36.19  by  Michael Scheer
*CMZ :  2.66/07 10/03/2010  09.50.10  by  Michael Scheer
*-- Author :    Michael Scheer   24/02/2010
      subroutine bunch(dtphase)
*KEEP,gplhint.
*KEND.

      use bunchmod

      implicit none

      double precision dtphase,weight,wlen,bunchar,bunlen,wscale
      double precision tpuls,emod,psi,dpsidp,omega,p0,r56,bessel,
     &  freqr,dum

      integer ibuff,nbuffp
      parameter (nbuffp=1000)
      real ran(1),ranbuff(nbuffp)

      integer ical,nelec,iwarn,ifail,n,nharm,mode,ncoef,nbuff

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,halbach.
      include 'halbach.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,cohfun.
      include 'cohfun.cmn'
*KEEP,trackf90.
      include 'trackf90.cmn'
*KEND.

      external coherentfun
      real coherentfun

      integer nfunp
      parameter (nfunp=200)
      real*8 fint(nfunp),coef(nfunp),xfun(nfunp)

      data ical/0/
      data iwarn/0/

      if (ical.eq.0) then

        if (ibunch.ne.0) then
          bunlen=bunchlen
          bunchar=bunchcharge
          nharm=nbunchharm
          ncoefcf=nbunchharm
          nelec=neinbunch
          p0=bunchp0
          r56=bunchr56
          mode=iubunch
        else if (iampli.ne.0) then
          bunlen=ampbunchlen
          bunchar=ampbunchcharge
          nharm=nampbunchharm
          ncoefcf=nampbunchharm
          nelec=-iamprep
          p0=ampbunchp0
          r56=ampbunchr56
          mode=iampcoh
        else
            write(lungfo,*)
     &        '*** Error in BUNCH: IBUNCH nor IAMPLI .ne. 0'
            write(6,*)
     &        '*** Error in BUNCH: IBUNCH nor IAMPLI .ne. 0'
            stop '*** Program WAVE aborted ***'
        endif

        modecf=mode

        if (p0.le.0.0d0) p0=1.0d-30
        if (r56.le.0.0d0) r56=1.0d-30

        dtphase=(wtra2i+(1.d0/(dmygamma*dmygamma))
     &    *(xstop-xstart)/2.d0)/clight1
        freqr=2.d0*pi1/dtphase*hbarev1

        if (ampfreq.eq.-9999.d0) then
          ampfreq=freqr
        endif  !(imampli.gt.0.and.ampfreq.eq.-9999.d0)

        omega=ampfreq/hbarev1
        wlen=wtoe1/ampfreq*1.0d-9

        tpuls=bunlen/clight1

        if (mode.eq.2) then

          dpsidp=twopi1/wlen*r56/dmyenergy
          wscale=0.0d0

          do n=1,nharm
            emod=exp(-(n*espread*dmyenergy*dpsidp)**2/2.0d0)
            call util_bessel(n,n*p0*dpsidp,bessel,ifail)
            if (ifail.ne.0.and.iwarn.eq.0) then
              write(lungfo,*)'*** Warning: Error in BUNCH from call to UTIL_BESSEL ***'
              write(lungfo,*)'*** Warning: further warnings will be suppressed ***'
              write(6,*)'*** Warning: Error in BUNCH from call to UTIL_BESSEL ***'
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

        else if (mode.eq.3) then

          open(unit=21,file='wave_phasespace.dat',status='old')

        else if (mode.eq.4) then

          ccf(0)=1.0d0
          open(unit=99,file='fourier-bunch.dat',status='old')
11        call util_skip_comment(99)
          read(99,*,end=99)dum
          ncoef=ncoef+1
          if (ncoef.ge.nharm) goto 99
          if (ncoef.gt.ncoefcfp) then
            write(lungfo,*)
     &        '*** Error in BUNCH: Number of coefficients on file'
            write(lungfo,*)'fourier-bunch.dat exceeds',ncoef,' ***'
            write(6,*)
     &        '*** Error in BUNCH: Number of coefficients on file'
            write(6,*)'fourier-bunch.dat exceeds',ncoef,' ***'
            stop '*** Program WAVE aborted ***'
          endif
          goto 11
99        rewind(99)

          if (ncoef.lt.nharm) then
            write(lungfo,*)
     &        '*** Warning in BUNCH: Number of coefficients on file'
            write(lungfo,*)
     &        'fourier-bunch.dat less than NAMPBUNCHHARM or NBUNCHHARM ***'
            write(6,*)
     &        '*** Warning in BUNCH: Number of coefficients on file'
            write(6,*)
     &        'fourier-bunch.dat less than NAMPBUNCHHARM or NBUNCHHARM ***'
            nharm=ncoef
          endif

          do n=1,nharm
            call util_skip_comment(99)
            read(99,*)ccf(n)
          enddo
          close(99)

          ncoefcf=nharm

        endif !mode

          ! Besser wscale nicht verwenden, da der inkohaerente Untergrund sonst
          ! skaliert wird. Besser Untergrund separat berechen und nachtraeglich
          ! addieren
        if (iampincoh.eq.0) then
          wscale=1.0d0
        endif

        if (mode.ne.3) then


          cfbunlen=bunlen
          cfxlenfou=cfbunlen
          ncoefcf=nharm
          print*,'--- Warning in subroutine bunch: util_random_func_init'
          print*,'--- (replacement for funlpx of CERNLIB) not yet tested'
          print*,'--- in this context. Be careful!'
          call util_random_func_init(coherentfun,nfunp,xfun,fint,coef)

        endif !iubunch

        ical=1

        ibuff=0
        if (nelec.gt.nbuffp) then
          nbuff=nbuffp
        else
          nbuff=nelec
        endif

      endif !ical

      if (mode.eq.2) then
c siehe Salding et al, NIM A 539 (2005) 499-526

        call ranmar(ran,1)
        dtphase=ran(1)*tpuls

        psi=omega*dtphase

        if (psi.ne.0.0d0.and.p0.ne.0.0d0.and.dpsidp.ne.0.0d0) then
          weight=0.0d0
          do n=1,nharm
            emod=exp(-(n*espread*dmyenergy*dpsidp)**2/2.0d0)
            call util_bessel(n,n*p0*dpsidp,bessel,ifail)
            if (ifail.ne.0.and.iwarn.eq.0) then
              write(lungfo,*)'*** Warning: Error in BUNCH from call to UTIL_BESSEL ***'
              write(lungfo,*)'*** Warning: further warnings will be suppressed ***'
              write(6,*)'*** Warning: Error in BUNCH from call to UTIL_BESSEL ***'
              write(6,*)'*** Warning: further warnings will be suppressed ***'
              iwarn=1
            endif
            weight=weight-2.0d0*emod*bessel*cos(n*psi)*wscale
          enddo
          weight=(1.0d0+weight)
     &      *bunchar/echarge1/nelec/wscale**2
        else
          weight=bunchar/echarge1/nelec
        endif

        else if (mode.eq.3) then

          call util_skip_comment(21)
          read(21,*)egamma,bunchx,xelec,yelec,zelec,ypelec,zpelec
          dtphase=bunchx/clight1

        else if (mode.eq.4) then

          ibuff=ibuff+1
          if (ibuff.gt.nbuff) ibuff=1
          if (ibuff.eq.1) then
          call util_random_func(nfunp,xfun,fint,coef,nbuff,ranbuff)
          endif

          dtphase=ranbuff(ibuff)/clight1

      endif !mode

      if (weight.lt.0.0d0) weight=0.0d0

      return
      end
