*CMZ :          31/12/2023  14.41.16  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.32.18  by  Michael Scheer
*CMZ :  4.01/03 17/05/2023  10.57.05  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  13.32.32  by  Michael Scheer
*CMZ :  4.01/00 22/02/2023  14.57.49  by  Michael Scheer
*-- Author : Michael Scheer
*KEEP,GPLHINT.
*KEND.
      program urad_phase_main

      use omp_lib
      use uradphasemod

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      double precision, dimension(:), allocatable :: z,y
      double precision, dimension(:,:), allocatable :: s

      double precision :: banwid=0.001,xbeta=0.0d0,
     &  perlen,shift,ebeam,curr,step,perl,
     &  pincen(3),pinw,pinh,park,wlen1,gamma,
     &  ephmin,ephmax,beffv,beffh,pherror,stosum(4),
     &  alphah,alphav,espread,harm,b0eff,rhv,
     &  betah,betav,eps0h,eps0v,pinx,piny,pinz,
     &  disph,dispph,dispv,disppv,bunchcharge,bunchlen,efi(3),bfi(3),rn(3),
     &  emith,emitv,pinxprop,pinwprop,pinhprop

      real xran(1),rr(2),axr,axi,ayr,ayi,azr,azi

      integer :: idebug=0,noranone,i,
     &  npiny,npinz,nper,nepho,modeph,modepin,modesphere,nharm,iy,iz,iobs,
     &  mthreads,nelec,icohere,ihbunch,ipho,iobph,iel,modebunch,ifieldprop,
     &  modewave=0,isto,nlpoi=0,nobsvprop,npinyprop,npinzprop

      namelist/uradphasen/
     &  perlen,shift,nper,beffv,beffh,
     &  ebeam,curr,step,noranone,nobsvprop,npinyprop,npinzprop,
     &  pinx,piny,pinz,pinw,pinh,npiny,npinz,modepin,modesphere,nharm,harm,
     &  nepho,ephmin,ephmax,pherror,pinxprop,pinwprop,pinhprop,
     &  mthreads,nelec,icohere,ihbunch,modeph,modebunch,ifieldprop,
     &  betah,betav,alphah,alphav,emith,emitv,espread,
     &  disph,dispph,dispv,disppv,bunchcharge,bunchlen

      integer :: irnsize=64,irnseed(64),ifixseed
      namelist/seedn/irnseed,ifixseed

      integer :: luna,istat,kalloc=1

      open(newunit=luna,file='urad_phase.nam',status='old',iostat=istat)
      if (istat.ne.0) then
        stop "*** Error: Could not open urad_phase.nam"
      endif

      read(luna,uradphasen)
      read(luna,seedn)

      close(luna)

      if(nelec.eq.1.and.noranone.eq.0) then
        noranone=1
        print*
        print*,'*** Changed NORANONE=0 to NORANONE=1, since NELEC=1'
        print*
      endif

      pincen=[pinx,piny,pinz]
      eps0h=emith*1.0d-9
      eps0v=emitv*1.0d-9

      bunchlen=bunchlen/1.0d9 !nm->m

      !print*,"sigz, sizp:",sigz/1000.0d0,sigzp/1000.0d0
      !print*,"sigy, siyp:",sigy/1000.0d0,sigyp/1000.0d0

      if (ifixseed.ne.0) then
        ifixseed=1
        call util_random_set_seed(irnsize,irnseed)
      endif

      if (nharm.gt.0.and.harm.gt.0.0d0) then

        gamma=ebeam/emassg1
        wlen1=wtoe1/abs(harm/nharm)
        perl=perlen/1000.0d0
        park=2.0d0*(wlen1/(perl*1.0d9/2.0d0/gamma**2)-1.0d0)

        if (park.lt.0.0d0) then
          write(6,*)
     &      '*** Error in urad_phase_main:'
          write(6,*)
     &      'Inconsistent values of nharm, harm, and perlen'
          write(6,*)' '
          stop
        endif

        park=sqrt(park)
        b0eff=park/(echarge1*perl/(2.*pi1*emasskg1*clight1))

        if (beffh.eq.0.0d0.and.beffv.ne.0d0) then
          beffv=beffv/abs(beffv)*b0eff
        else if (beffv.eq.0.0d0.and.beffh.ne.0d0) then
          beffh=beffh/abs(beffh)*b0eff
        else
          rhv=beffh/beffv
          beffh=b0eff/sqrt(1.0d0+1.0d0/rhv**2)*beffh/abs(beffh)
          beffv=beffh/rhv
        endif

      endif

      npiny=max(1,npiny)
      npinz=max(1,npinz)

      open(newunit=luna,file='urad_phase.pin')
      write(luna,*)npinz,npiny,pinw,pinh
      write(luna,*)pincen
      close(luna)

      if (mthreads.lt.0) then
        mthreads=OMP_GET_MAX_THREADS()
      else if (mthreads.eq.0) then
        mthreads=1
      endif

      call urad_phase(
     &  mthreads,nelec,noranone,icohere,modebunch,bunchlen,bunchcharge,ihbunch,
     &  perlen,shift,nper,beffv,beffh,
     &  ebeam,curr,step,nlpoi,
     &  pincen,pinw,pinh,npiny,npinz,modepin,modesphere,
     &  nepho,ephmin,ephmax,banwid,
     &  xbeta,betah,alphah,betav,alphav,espread,emith,emitv,
     &  disph,dispph,dispv,disppv,
     &  modeph,pherror,modewave
     &  )

      open(newunit=luna,file='urad_phase.fld')

      do iobs=1,nobsv_u
        do ipho=1,nepho_u
          iobph=iobs+nobsv_u*(ipho-1)

          rn(1)=real(arad_u(2,iobph)*conjg(arad_u(6,iobph))-arad_u(3,iobph)*conjg(arad_u(5,iobph)))
          rn(2)=real(arad_u(3,iobph)*conjg(arad_u(4,iobph))-arad_u(1,iobph)*conjg(arad_u(6,iobph)))
          rn(3)=real(arad_u(1,iobph)*conjg(arad_u(5,iobph))-arad_u(2,iobph)*conjg(arad_u(4,iobph)))
          rn=rn/norm2(rn)

          axr=real(arad_u(1,iobph))
          axi=imag(arad_u(1,iobph))
          ayr=real(arad_u(2,iobph))
          ayi=imag(arad_u(2,iobph))
          azr=real(arad_u(3,iobph))
          azi=imag(arad_u(3,iobph))

          write(luna,'(3(1pe15.6e3),i10,21(1pe15.6e3))')
     &      obsv_u(1:3,iobs),ipho,epho_u(ipho),stokes_u(1:4,iobph),pow_u(iobs),
     &      real(arad_u(1,iobph)),imag(arad_u(1,iobph)),
     &      real(arad_u(2,iobph)),imag(arad_u(2,iobph)),
     &      real(arad_u(3,iobph)),imag(arad_u(3,iobph)),
     &      real(arad_u(4,iobph)),imag(arad_u(4,iobph)),
     &      real(arad_u(5,iobph)),imag(arad_u(5,iobph)),
     &      real(arad_u(6,iobph)),imag(arad_u(6,iobph)),
     &      rn

        enddo
      enddo
      close(luna)

      allocate(z(nobsv_u),y(nobsv_u),s(npinz_u,npiny_u))

      open(newunit=luna,file='urad_phase.flx')

      z=obsv_u(3,1:npinz_u)

      do iy=1,npiny_u
        iobs=npinz_u*(iy-1)+iy
        y(iy)=obsv_u(2,iobs)
      enddo

      do ipho=1,nepho_u
        if (modepin.eq.0.and.npinz_u.ge.3.and.npiny_u.ge.3) then
          do isto=1,4
            iobs=0
            do iz=1,npinz_u
              do iy=1,npiny_u
                iobs=iobs+1
                iobph=iobs+nobsv_u*(ipho-1)
                s(iz,iy)=stokes_u(isto,iobph)
              enddo
            enddo
            call util_spline_integral_2d(npinz_u,npiny_u,z,y,s,stosum(isto),
     &        istat,kalloc)
            kalloc=0
          enddo !isto
          write(luna,*)ipho,epho_u(ipho),stosum
        else
          iobs=0
          stosum=0.0d0
          do iz=1,npinz_u
            do iy=1,npiny_u
              iobs=iobs+1
              iobph=iobs+nobsv_u*(ipho-1)
              stosum=stosum+stokes_u(1:4,iobph)
            enddo
          enddo
          write(luna,*)ipho,epho_u(ipho),stosum/nobsv_u*pinw*pinh
        endif
      enddo
      close(luna)

      open(newunit=luna,file='urad_phase.bun')
      if (ihbunch.le.0) then
        write(luna,*)
      else
        do iel=1,nelec_u/ihbunch_u*nepho_u
          if(fbunch_u(21,iel).ne.0.0d0) then
            write(luna,*)fbunch_u(:,iel)
          endif
        enddo
      endif
      close(luna)

      call  util_random_get_seed(irnsize,irnseed)

      open(newunit=luna,file='urad_phase.seeds',status='unknown')
      write(luna,*)irnsize
      do i=1,irnsize
        write(luna,*)i,irnseed(i)
      enddo
      flush(luna)
      close(luna)

      if (ifieldprop.ne.0) then
        pinxprop_u=pinxprop
        pinwprop_u=pinwprop
        pinhprop_u=pinhprop
        npinyprop_u=max(1,npinyprop)
        npinzprop_u=max(1,npinzprop)

        call urad_phase_prop(mthreads)

        open(newunit=luna,file='urad_phase.fdp')

        do iobs=1,nobsvprop_u
          do ipho=1,nepho_u
            iobph=iobs+nobsvprop_u*(ipho-1)

            rn(1)=real(aradprop_u(2,iobph)*conjg(aradprop_u(6,iobph))-aradprop_u(3,iobph)*conjg(aradprop_u(5,iobph)))
            rn(2)=real(aradprop_u(3,iobph)*conjg(aradprop_u(4,iobph))-aradprop_u(1,iobph)*conjg(aradprop_u(6,iobph)))
            rn(3)=real(aradprop_u(1,iobph)*conjg(aradprop_u(5,iobph))-aradprop_u(2,iobph)*conjg(aradprop_u(4,iobph)))
            rn=rn/norm2(rn)

            axr=real(aradprop_u(1,iobph))
            axi=imag(aradprop_u(1,iobph))
            ayr=real(aradprop_u(2,iobph))
            ayi=imag(aradprop_u(2,iobph))
            azr=real(aradprop_u(3,iobph))
            azi=imag(aradprop_u(3,iobph))

            write(luna,'(3(1pe15.6e3),i10,20(1pe15.6e3))')
     &        obsvprop_u(1:3,iobs),ipho,epho_u(ipho),stokesprop_u(1:4,iobph),
     &        real(aradprop_u(1,iobph)),imag(aradprop_u(1,iobph)),
     &        real(aradprop_u(2,iobph)),imag(aradprop_u(2,iobph)),
     &        real(aradprop_u(3,iobph)),imag(aradprop_u(3,iobph)),
     &        real(aradprop_u(4,iobph)),imag(aradprop_u(4,iobph)),
     &        real(aradprop_u(5,iobph)),imag(aradprop_u(5,iobph)),
     &        real(aradprop_u(6,iobph)),imag(aradprop_u(6,iobph)),
     &        rn

          enddo
        enddo

        close(luna)

      endif

      end
