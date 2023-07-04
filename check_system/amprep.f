*CMZ :  4.00/15 31/05/2022  09.05.33  by  Michael Scheer
*CMZ :  4.00/11 28/06/2021  10.33.06  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine amprep

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.
      use sourcef90
      use observf90
      use afreqf90
      use bunchmod

      implicit none

      double precision, dimension (:), allocatable :: wsspec
      double precision, dimension (:,:), allocatable :: wsstokes

      real, dimension (:), allocatable :: pherr
      real, dimension(:), allocatable :: phiran

      real xran(5),pran(2),rr(2)

      double complex apol,amp0(3),damp(3),amp(3),zexp,
     &  apolh,apolr,apoll,apol45,phbu

      double precision :: dtelec0,dtelec,dtpho,t0,perlen,dph,dobs(3),drn(3),
     &  cosang,ang,
     &  dobsn(3),r0(3),r(3),dr(3),t,dt,dist,obs(3),dlam,om,om1,fd,dist0,dl,dobsbu(3),
     &  stok1,stok2,stok3,stok4,speck,sqnbunch,sqnphsp,pow,perr,obs0(3),
     &  x2,y2,z2,vx2,vy2,vz2,dgamma,eix,eiy,eiz,efx,efy,efz,vn,vf0,v0,v,obscen(3),
     &  xi,yi,zi,vxi,vyi,vzi,xe,ye,ze,vxe,vye,vze,bshift=0.5d0,gamma,
     &  alphah,alphav,beta0h,beta0v,dpp,zpp,zz,zpi,ypi,zzp,zp,yp,yyp,yy,
     &  sigh,sigph,sigv,sigpv,s0v,s0h,gammav,gammah,fillb(29),zpmax,
     &  deflpar,xkx,pathlen,tof,dtdpp,t1,t2,distcorr,b0par(3),xb0par(3),
     &  b0opt,xb0opt,ab0(3),b0p(3),ds,upow,udgamtot,flow,fhigh,dum(3),dist00,
     &  obangh0,obangh,obangv0,obangv,drn00(3),drn0(3),dr0(3),r00(3),dr00(3),
     &  rm(3),dobsm(3),dx,beta,rpin,ppin,beff,park,wlen1,eharm1

      double precision dppv(13),dtev(13),dpp2(13),ws1(13),ws2(13),ws3(13),ws4(13)

      double complex, dimension (:), allocatable :: uampx,uampy,uampz
      double precision, dimension (:,:), allocatable :: utraxyz,ustokes

      integer :: kfreq,iobsv,nper,i,np2,nelec,mbunch,meinbunch,ibu,kran,ib0max,k,
     &  ifail,ndimu,nstepu,l,ith

      integer idatetime(8)
      character(10) dtday,dttime,dtzone
      character(3) cmodph

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,track0.
      include 'track0.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,observ.
      include 'observ.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,optic.
      include 'optic.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,ellip.
      include 'ellip.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.
c+seq,uservar.

      stop "*** amprep is obsolete: Use amprep_omp etc. ***"
      call date_and_time(dtday,dttime,dtzone,idatetime)

      write(cmodph,'(I3)') modeph

      write(6,*)
c      write(6,*)'     Starting calculations in AMPREP with  MODEPH =' // cmodph //': '
      write(6,*)'     Starting calculations in AMPREP: '
     &  ,dttime(1:2),':',dttime(3:4),':',dttime(5:6)
      write(6,*)

      r00=[x0,y0,z0]
      dr00=[xf0-x0,yf0-y0,zf0-z0]
      r00=r00+dr00/2.0d0

      perlen=norm2(dr00)
      drn00=dr00/perlen

      vxi=vx0
      vyi=vy0
      vzi=vz0

      r0=r00
      drn0=drn00
      r=r0
      drn=drn0

      v0=norm2([vx0,vy0,vz0])
      efx=vx0/v0
      efy=vy0/v0
      efz=vz0/v0

      vf0=norm2([vxf0,vyf0,vzf0])
      efx=vxf0/vf0
      efy=vyf0/vf0
      efz=vzf0/vf0

      ds=dtim0*vf0
      ndimu=nco*1.1

      allocate(uampx(nfreq),uampy(nfreq),uampz(nfreq),
     &  utraxyz(14,ndimu),ustokes(4,nfreq))

      if (nfreq.eq.1) then
        flow=freq(1)
        fhigh=freq(1)
      else
        flow=freqlow
        fhigh=freqhig
      endif

      beff=sqrt(bymx**2+bzmx**2)
      park=echarge1*beff*perlen/(2.*pi1*emasskg1*clight1)
      wlen1=(1+park**2/2.)/2./dmygamma**2*perlen*1.0d9

      if (wlen1.ne.0.0) then
        eharm1=wtoe1/wlen1
      else
        eharm1=0.0d0
      endif

      dtpho=perlen/clight1

      nper=iabs(kampli)

      allocate(pherr(nper),affe(3,nfreq*nobsv),phiran(max(1,nbunch)),
     &  wsspec(nfreq*nobsv),wsstokes(4,nfreq*nobsv))

      np2=nper/2

      affe=(0.0D0,0.0D0)
      spec=0.0d0
      stokes=0.0d0

      call util_random_gauss(nper,pherr,rr)

      if (eharm1.ne.0.0d0) then
        om1=eharm1/hbarev1 ! omega = 2 pi nu = 2 pi / T; 2 pi / omega = T
        pherr=pherr*pherror/360.0d0*twopi1/om1
      else
        pherr=0.0d0
      endif

      mbunch=max(1,nbunch)
      meinbunch=max(1,neinbunch)
      nelec=max(1,mbunch*meinbunch)

      if (ibunch.ne.0.and.meinbunch.ne.1) then
        stop "*** Error in amprep: Neinbunch not one is not yet implemented ***"
      endif

      if (ibunch.ne.0.and.bunchcharge.ne.0.0d0) then
        sqnbunch=mbunch
        sqnphsp=sqrt(bunchcharge/echarge1)
     &    *meinbunch
     &    /(bunchcharge/echarge1)
        bunnor=1.0d0/mbunch
      else
        sqnbunch=mbunch
        sqnphsp=sqrt(dble(meinbunch))
        bunnor=1.0d0/mbunch
      endif

      call util_random(nbunch,phiran)
      phiran=phiran*twopi1
      if (ibunch.eq.-1) phiran(1)=0.0

      if (ibunch.ne.0) then
        if (iubunch.eq.1) then
          alphah=-betaph/2.0d0
          gammah=(1.0d0+alphah**2)/betah
          beta0h=1.0d0/gammah
          s0h=alphah/gammah
          alphav=-betapv/2.0d0
          gammav=(1.0d0+alphav**2)/betav
          beta0v=1.0d0/gammav
          s0v=alphav/gammah
          sigh=sqrt(eps0h*beta0h)
          sigph=sqrt(eps0h/beta0h)
          sigv=sqrt(eps0v*beta0v)
          sigpv=sqrt(eps0v/beta0v)
        else if (iubunch.ne.-1.and.iubunch.ne.0.and.
     &      (noemitph.eq.0.or.noespreadph.ne.0)) then
          stop "*** Error in AMPREP: IUBUNCH must be one of [0,1,-1] ***"
        endif
      endif

      ielec=0

      do ibu=1,nelec

        isub=1 ! später Schleife über neinbunch

        wsspec=0.0d0
        wsstokes=0.0d0

        ielec=ielec+1
        phbu=cdexp(dcmplx(0.0d0,phiran(ibu)))
        bunchx=0.0d0

        xi=x0
        yi=y0
        zi=z0

        zpi=vz0/vx0
        ypi=vy0/vx0

        x2=xf0
        y2=yf0
        z2=zf0

        vx2=vxf0
        vy2=vyf0
        vz2=vzf0

        dpp=0.0d0
        gamma=dmygamma
        dtelec=dtelec0

        if (ibunch.ne.0.and.(ibunch.ne.-1.or.ielec.gt.1)) then

          kran=0
          if (noespreadph.eq.0) kran=1
          if (noemitph.eq.0) kran=5

          call util_random_gauss(kran,xran,rr)

          if (iamppin.eq.3) then
            call util_random(2,pran)
            obs(1)=pincen(1)
            if (iamppincirc.eq.0) then
              obs(2)=pincen(2)+(pran(1)-0.5)*pinw
              obs(3)=pincen(3)+(pran(2)-0.5)*pinh
            else
              rpin=(pran(1)-0.5)*pinr
              ppin=pran(2)*twopi1
              obs(2)=pincen(2)+rpin*cos(ppin)
              obs(3)=pincen(3)+rpin*sin(ppin)
            endif
          endif

          if (noespreadph.eq.0) then
            dpp=espread*xran(1)
            gamma=(1.0d0+dpp)*dmygamma
          endif

          ! assume beta(s)=beta0(s)+s**2/beta(0) and alpha0=-s/beta(0)
          ! and a drift transfer-matrix ((1,s),(1,0))

          if (noemitph.eq.0) then

            if (iubunch.eq.0) then

              zz=bsigz(1)*xran(2)
              zzp=bsigzp(1)*xran(3)

              zi=zz-x0*zzp !inverse transformation
              zpi=zzp

              yy=bsigy(1)*xran(4)
              yyp=bsigyp(1)*xran(5)

              yi=yy-x0*yyp
              ypi=yyp

            else if (iubunch.eq.1) then

              zz=sigh*xran(2)
              zzp=sigph*xran(3)

              zi=zz-s0h*zzp !inverse transformation
              zpi=zzp

              yy=sigv*xran(4)
              yyp=sigpv*xran(5)

              yi=yy-s0v*yyp
              ypi=yyp

            endif !(noemitph.eq.0)

          else if (iubunch.eq.-1) then

            call ubunch(xi,yi,zi,ypi,zpi,gamma,dt)
            dpp=(gamma-dmygamma)/dmygamma
            gamma=(1.0d0+dpp)*dmygamma

          endif !iubunch

          if (noemitph.eq.0.or.noespreadph.eq.0) then

            ! simple treatment of closed orbit, assume small angles

            zi=zi+z0
            zpi=zpi+vz0/vx0

            yi=yi+y0
            ypi=ypi+vy0/vx0

            vn=clight1*dsqrt((1.0d0-1.0d0/gamma)*(1.0d0+1.0d0/gamma))

            zi=zi+dpp*disp0
            zpi=zpi+dpp*ddisp0

            vxi=vn/sqrt(1.0d0+ypi**2+zpi**2)
            vyi=vxi*ypi
            vzi=vxi*zpi

          endif !(noemitph.eq.0.or.noespreadph.eq.0)

        endif !(ibunch.ne.-1.or.ielec.ne.1)

        vn=norm2([vxi,vyi,vzi])
        eix=vxi/vn
        eiy=vyi/vn
        eiz=vzi/vn

        x2=x0+eix*xlellip
        y2=y0+eiy*xlellip
        z2=z0+eiz*xlellip

        r0=[xi,yi,zi]
        dr0=[x2-xi,y2-yi,z2-zi]
        r0=r0+dr0/2.0d0

        do iobsv=1,nobsv

          if (iamppin.ne.3) then
            obs=obsv(1:3,iobsv)
          endif

          call urad_omp(icharge,
     &      gamma,udgamtot,
     &      xi,yi,zi,vxi,vyi,vzi,
     &      xf0,yf0,zf0,efx,efy,efz,
c     &      x0+eix*xlellip,y0+eiy*xlellip,z0+eiz*xlellip,eix,eiy,eiz,
     &      x2,y2,z2,vx2,vy2,vz2,dtelec,ds,
     &      0,nstepu,ndimu,utraxyz,
     &      obs(1),obs(2),obs(3),flow,fhigh,
     &      nfreq,freq,uampx,uampy,uampz,ustokes,upow,
     &      ieneloss,ivelofield,ifail)

          r0=[xi,yi,zi]
          dr0=[x2-xi,y2-yi,z2-zi]
          r0=r0+dr0/2.0d0

          if (ielec.eq.1) then
            pow=specpow(iobsv)
            specpow(iobsv)=0.0d0
          endif

          do kfreq=1,nfreq

            ifrob=kfreq+nfreq*(iobsv-1)
            iobfr=iobsv+nobsv*(kfreq-1)

            om=freq(kfreq)/hbarev1

            amp0=[uampx(kfreq),uampy(kfreq),uampz(kfreq)]*1.0d3/sqrt(specnor/dmycur*0.10d0) !urad
            amp=(0.0d0,0.0d0)
            t=0.0d0

            do i=1,nper

              r=r0+(i-np2-1)*dr0
              dobs=obs-r
              dist0=norm2(obs-r0)
              dist=norm2(dobs)

              dt=xlellip/clight1*((1.0d0+parkell**2/2.0d0)/2.0d0/gamma**2+
     &          (((ypi-dobs(2)/dobs(1))**2+(zpi-dobs(3)/dobs(1))**2))/2.0d0)

c              amp0(1)=complex(reaima(1,1,iobfr),reaima(1,2,iobfr))
c              amp0(2)=complex(reaima(2,1,iobfr),reaima(2,2,iobfr))
c              amp0(3)=complex(reaima(3,1,iobfr),reaima(3,2,iobfr))

              t=t+dt
              dph=om*(t+pherr(i))

              zexp=cdexp(dcmplx(0.0d0,dph))
              damp=amp0*zexp*dist0/dist
              amp=amp+damp

              if (ihbunch.ne.0) then

                if (iobsv.eq.icbrill.and.(ibu.eq.1.or.mod(ielec,ihbunch).eq.0)) then

                  if (i.eq.1) then
                    fillb(5)=r(1)
                    fillb(6)=r(2)
                    fillb(7)=r(3)
                    fillb(8)=ypi
                    fillb(9)=zpi
                  else if (i.eq.nper) then
                    fillb(10:12)=r+dr
                    fillb(13)=ypi
                    fillb(14)=zpi
                  endif

                endif
              endif

              if (kfreq.eq.1.and.ibu.eq.1)
     &          specpow(iobsv)=specpow(iobsv)+pow*(dist0/dist)**2

            enddo !nper

            if(speccut.gt.0.0d0) then
              if(freq(kfreq).gt.speccut*ecdipev1*dmyenergy**2*ecmax(1)) then
                amp=(0.0d0,0.0d0)
              endif
            endif

            amp=amp*reflec

            if (ipola.eq.0) then
              speck=sum(amp*conjg(amp))*specnor*bunnor
            else    !ipola
              apol=
     &           amp(1)*conjg(vpola(1))
     &          +amp(2)*conjg(vpola(2))
     &          +amp(3)*conjg(vpola(3))
              speck=dreal(apol*conjg(apol))*specnor*bunnor
            endif   !ipola

            wsspec(iobfr)=wsspec(iobfr)+speck

            if (istokes.ne.0) then

              apolh=
     &          amp(1)*conjg(vstokes(1,1))
     &          +amp(2)*conjg(vstokes(1,2))
     &          +amp(3)*conjg(vstokes(1,3))

              apolr=
     &          amp(1)*conjg(vstokes(2,1))
     &          +amp(2)*conjg(vstokes(2,2))
     &          +amp(3)*conjg(vstokes(2,3))

              apoll=
     &          amp(1)*conjg(vstokes(3,1))
     &          +amp(2)*conjg(vstokes(3,2))
     &          +amp(3)*conjg(vstokes(3,3))

              apol45=
     &          amp(1)*conjg(vstokes(4,1))
     &          +amp(2)*conjg(vstokes(4,2))
     &          +amp(3)*conjg(vstokes(4,3))

              stok1=apolr*conjg(apolr)+apoll*conjg(apoll)
              stok2=-stok1+2.0d0*apolh*conjg(apolh)
              stok3=2.0d0*apol45*conjg(apol45)-stok1
              stok4=apolr*conjg(apolr)-apoll*conjg(apoll)

              wsstokes(1,iobfr)=wsstokes(1,iobfr)+stok1*specnor*bunnor
              wsstokes(2,iobfr)=wsstokes(2,iobfr)+stok2*specnor*bunnor
              wsstokes(3,iobfr)=wsstokes(3,iobfr)+stok3*specnor*bunnor
              wsstokes(4,iobfr)=wsstokes(4,iobfr)+stok4*specnor*bunnor

            endif !istokes

            spec(iobfr)=spec(iobfr)+speck
            stokes(1:4,iobfr)=stokes(1:4,iobfr)+wsstokes(1:4,iobfr)
            affe(:,ifrob)=affe(:,ifrob)+phbu*amp

            if (ihbunch.ne.0) then
              if (iobsv.eq.icbrill.and.(ibu.eq.1.or.mod(ielec,ihbunch).eq.0)) then

                fillb(1)=ibun
                fillb(2)=isub
                fillb(3)=ielec
                fillb(4)=bunchx
                fillb(15)=gamma*emassg1
                fillb(16)=(gamma+dgamma)*emassg1
                fillb(17)=obsv(1,iobsv)
                fillb(18)=obsv(2,iobsv)
                fillb(19)=obsv(3,iobsv)
                fillb(20)=kfreq
                fillb(21)=freq(kfreq)
                fillb(22)=wsspec(iobfr)*nelec

                if (istokes.ne.0) then
                  fillb(23)=wsstokes(1,iobfr)*nelec
                  fillb(24)=wsstokes(2,iobfr)*nelec
                  fillb(25)=wsstokes(3,iobfr)*nelec
                  fillb(26)=wsstokes(4,iobfr)*nelec
                else
                  fillb(23)=fillb(22)
                  fillb(24:26)=0.0d0
                endif !istokes

                fillb(27)=specpow(iobsv)
                fillb(28)=1
                fillb(29)=dtelec

                call hfm(nidbunch,fillb)

              endif
            endif

          enddo !kfreq

        enddo !nobsv

      enddo !nbunch

      do iobsv=1,nobsv
        do kfreq=1,nfreq
          ifrob=kfreq+nfreq*(iobsv-1)
          iobfr=iobsv+nobsv*(kfreq-1)
          reaima(1:3,1,iobfr)=dreal(affe(1:3,ifrob))/sqnbunch
          reaima(1:3,2,iobfr)=dimag(affe(1:3,ifrob))/sqnbunch
        enddo !kfreq
      enddo !nobsv

      call date_and_time(dtday,dttime,dtzone,idatetime)

      write(6,*)
c      write(6,*)'     Finishing calculations in AMPREP with  MODEPH =' // cmodph //': '
      write(6,*)'     Finishing calculations in AMPREP: '
     &  ,dttime(1:2),':',dttime(3:4),':',dttime(5:6)
      write(6,*)

      deallocate(uampx,uampy,uampz,utraxyz,ustokes)
      deallocate(pherr,affe,phiran,wsspec,wsstokes)

      return
      end
