*CMZ :          13/08/2024  10.11.51  by  Michael Scheer
*CMZ :  4.01/05 26/04/2024  10.49.56  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase(
     &  mthreads,nelec,noranone,icohere,modebunch,bunchlen,bunchcharge,ihbunch,
     &  perlen,shift,nper,beffv,beffh,
     &  ebeam,curr,step,nlpoi,
     &  pincen,pinw,pinh,npiny,npinz,modepin,modesphere,
     &  nepho,ephmin,ephmax,banwid,
     &  xbeta,betah,alphah,betav,alphav,espread,emith,emitv,
     &  disph,dispph,dispv,disppv,
     &  modeph,pherror,phgshift,modewave
     &  )

      use omp_lib
      use uradphasemod

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.
c+seq,uservar.

      double complex :: rea(3),expsh

      double precision
     &  perlen,shift,ebeam,curr,step,banwid,
     &  pincen(3),pinw,pinh,betah,alphah,betav,alphav,
     &  ephmin,ephmax,beffv,beffh,pherror,phgshift,espread,emith,emitv,
     &  disph,dispph,dispv,disppv,y,z,dy,dz,ymin,zmin,bunchlen,bunchcharge,
     &  xbeta,df,xx,yy,zz,r,xn,yn,zn,h2

      integer
     &  npiny,npinz,nper,nepho,mthreads,nelec,icohere,ihbunch,i,nlpoi,
     &  modeph,modepin,modesphere,modebunch,iy,iz,iobsv,noranone,modewave,
     &  icbrill,iobs,iobfr,ifrq

c      if (modewave.ne.0) call util_zeit_kommentar(6,'Entered urad_phase')
      call util_zeit_kommentar_delta(6,'Entered urad_phase',1)

      mthreads_u=mthreads

      nelec_u=nelec
      noranone_u=noranone
      icohere_u=icohere
      modebunch=modebunch_u
      bunchlen_u=bunchlen
      bunchcharge_u=bunchcharge
      ihbunch_u=ihbunch

      perlen_u=perlen/1000.0d0
      shift_u=shift/1000.0d0
      nper_u=nper
      beffv_u=beffv
      beffh_u=beffh

      ebeam_u=ebeam
      gamma_u=ebeam_u/emassg1
      step_u=step/1000.0d0
      nstep_u=max(1,nint(perlen_u/step_u))

      curr_u=curr

      pincen_u=pincen/1000.0d0
      pinw_u=pinw/1000.0d0
      pinh_u=pinh/1000.0d0
      npiny_u=npiny
      npinz_u=npinz
      modepin_u=modepin

      ephmin_u=ephmin_u
      ephmax_u=ephmax_u
      banwid_u=banwid
      nepho_u=nepho

      npiny_u=max(1,npiny_u)
      npinz_u=max(1,npinz_u)

      nlpoi_u=nlpoi

      if (modepin.ne.1) then
        nobsv_u=npiny_u*npinz_u
      else
        npinz_u=1
        npiny_u=1
        nobsv_u=1
      endif

      allocate(epho_u(nepho),obsv_u(3,nobsv_u),
     &  arad_u(6,nobsv_u*nepho_u),
     &  specpow_u(nobsv_u),
     &  stokes_u(4,nobsv_u*nepho_u),pow_u(nobsv_u)
     &  )

      if (ihbunch_u.gt.0) then
        allocate(fbunch_u(41,nelec_u/ihbunch_u*nepho_u))
        fbunch_u=0.0d0
      else if (ihbunch_u.lt.0) then
        allocate(fbunch_u(41,nobsv_u*nelec_u/(-ihbunch_u)*nepho_u))
        fbunch_u=0.0d0
      endif

      stokes_u=0.0d0
      specpow_u=0.0d0
      arad_u=(0.0d0,0.0d0)
      pow_u=0.0d0

      if (npiny_u.eq.1) then
        dy=0.0d0
        ymin=pincen_u(2)
      else
        dy=pinh_u/(npiny_u-1)
        ymin=pincen_u(2)-pinh_u/2.0d0
      endif

      if (npinz_u.eq.1) then
        dz=0.0d0
        zmin=pincen_u(3)
      else
        dz=pinw_u/(npinz_u-1)
        zmin=pincen_u(3)-pinw_u/2.0d0
      endif

      iobsv=0
      y=ymin-dy
      do iy=1,npiny_u
        y=y+dy
        z=zmin-dz
        do iz=1,npinz_u
          iobsv=iobsv+1
          z=z+dz
          obsv_u(1,iobsv)=pincen_u(1)
          obsv_u(2,iobsv)=y
          obsv_u(3,iobsv)=z
          if (modesphere.ne.0) then
            !all util_break
            xx=obsv_u(1,iobsv)
            yy=obsv_u(2,iobsv)
            zz=obsv_u(3,iobsv)
c            r=sqrt(xx*xx+yy*yy+zz*zz)
            h2=(zz**2+yy**2)/xx**2
            if (h2.lt.0.01) then
c              r=xx*(1.0d0+h2/2.0d0-h2**2/8.0d0)
              r=xx*(1.0d0+(((((-0.0205078125D0*h2+0.02734375D0)*h2
     &      -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2)
            else
              r=xx*(1.0d0+sqrt(1.0d0+h2))
            endif
            xn=xx/r
            yn=yy/r
            zn=zz/r
            obsv_u(1,iobsv)=xn*xx
            obsv_u(2,iobsv)=yn*xx
            obsv_u(3,iobsv)=zn*xx
          endif
        enddo
      enddo

      nepho_u=max(1,nepho_u)
      if (nepho_u.gt.1) then
        df=(ephmax-ephmin)/(nepho_u-1)
        do i=1,nepho_u
          epho_u(i)=ephmin+(i-1)*df
        enddo
      else
        epho_u(1)=(ephmin+ephmax)/2.0d0
      endif

      xbeta_u=xbeta
      betah_u=betah
      alphah_u=alphah
      betav_u=betav
      alphav_u=alphav
      emith_u=emith
      emitv_u=emitv
      disph_u=disph
      dispph_u=dispph
      dispv_u=dispv
      dispph_u=disppv
      espread_u=espread

      modeph=modeph_u
      pherror_u=pherror
      phgshift_u=phgshift

c      call urad_spline(modewave)
c      stop
c      if (modewave.eq.2) then
c        call urad_nnb(modewave)
c      else if (modewave.eq.3) then
c        call urad_spline(modewave)
c      else
      call urad_amprep(modewave)
c      endif

      stokes_u=stokes_u/1.0d6 ! photons/mm**2
      if (ihbunch.ne.0) then
        fbunch_u(4:14,:)=fbunch_u(4:14,:)*1000.0d0 ! mm
        fbunch_u(17:19,:)=fbunch_u(17:19,:)*1000.0d0 ! mm
        fbunch_u(22:26,:)=fbunch_u(22:26,:)/1.0d6 ! 1/mm**2
      endif

      icbrill=nobsv_u/2+1
      if (abs(phgshift).eq.9999.0d0) then
        do ifrq=1,nepho
          iobfr=icbrill+nobsv_u*(ifrq-1)
          rea(1:2)=(0.0d0,0.0d0)
          rea(3)=arad_u(3,iobfr)
          expsh=rea(3)/abs(rea(3))*1.0d3
          if (phgshift.eq.-9999.0d0) expsh=expsh*cdexp(dcmplx(0.0d0,-pi1/2.0d0))
          DO iobs=1,nobsv_u
            iobfr=iobs+nobsv_u*(ifrq-1)
            arad_u(1:6,iobfr)=arad_u(1:6,iobfr)/expsh
          enddo
        enddo
      else if (phgshift.ne.0.0d0) then
        expsh=cdexp(dcmplx(0.0d0,phgshift))*1.0d3
        arad_u=arad_u/expsh
      else
        arad_u=arad_u/1.0d3
      endif

      obsv_u=obsv_u*1000.0d0

c      if (modewave.ne.0) call util_zeit_kommentar(6,'Leaving urad_phase')
      call util_zeit_kommentar_delta(6,'Leaving urad_phase',0)

      end
