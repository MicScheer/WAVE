*CMZ :          21/08/2026  16.04.50  by  Michael Scheer
*-- Author :    Michael Scheer   17/08/2026
      SUBROUTINE GENPHO
*KEEP,GPLHINT.
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

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,phasef90u.
      include 'phasef90u.cmn'
*KEEP,phasewsf90u.
      include 'phasewsf90u.cmn'
*KEEP,wbetaf90u.
      include 'wbetaf90u.cmn'
*KEEP,wbetaf90u.
      include 'wbetaf90u.cmn'
*KEND.

      use ompmod
      use omp_lib
      use wobsvmod

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,reargf90.
      include 'reargf90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,wbetaf90.
      include 'wbetaf90.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,genpho.
      include 'genpho.cmn'
*KEND.

      complex*16, dimension(:,:), allocatable :: arad
      complex*16 efc(3),bfc(3),expsh,rea(8)

      real*8 :: specnor_si,
     &  zi=0.0d0, !dummy
     &  yi=0.0d0 !dummy

      real*8
     &  s0h,beta0h,gamma0h,
     &  s1,s2,beta1h,betap1h,alpha1h,gamma1h,phase1h,
     &  s2h,beta2h,betap2h,alpha2h,gamma2h,phase2h,
     &  s0v,beta0v,gamma0v,
     &  s1v,beta1v,betap1v,alpva1v,gamma1v,phase1v,
     &  s2v,beta2v,betap2v,alpva2v,gamma2v,phase2v,obs(3,mobsv)

      real, dimension(:), allocatable :: photons,electrons

      real :: ephomin,ephomax,
     &  sigz,sigzp,sigy,sigyp,z,y,g=1.0,zmin,zmax,ymin,ymax

      integer :: iz,iy,nedz,nedy,ifrq,iobs,ngam,iel,l,i,lunapho,lunaele,iepho,ndimpho,
     &  jobs,jobfr,kpin,ndimele,
     &  ndimapho=9,ndimaele=4

      character(12) chspacer

      if (moderan.gt.0) then
        ngam=npho
      else
        ngam=mobsv
      endif

      ndimpho=ndimapho*ngam*nfreq*NELECGENPHO
      ndimele=ndimaele*nelecgenpho
      print*,"GENPHO: ndimpho:",ndimpho
      allocate(photons(ndimpho),electrons(ndimele))
      photons=0.0
      electrons=0.0
      allocate(arad(6,mobsv*nfreq))
      arad=(0.0d0,0.0d0)

      s1=sourcea(1,1,1)
      s2=sourcee(1,1,1)

      call util_beta_function_drift(
     &  s0h,beta0h,gamma0h,
     &  s1,betah,betaph,alpha1h,gamma1h,phase1h,
     &  s2h,beta2h,betap2h,alpha2h,gamma2h,phase2h)

      call util_beta_function_drift(
     &  s0v,beta0v,gamma0v,
     &  s1,betav,betapv,alpva1v,gamma1v,phase1v,
     &  s2v,beta2v,betap2v,alpva2v,gamma2v,phase2v)

c      print*,''
c      print*,'      --- Subroutine GENPHO ---   '
c      print*,''

      if (
     &    abs((s2+s1)/2.0d0).gt.1.0d0/dble(myinum)
     &    ) then
        print*,''
        print*,"*** Warning in GENPHO: X of source center not in origin ***"
        print*,''
      endif

      if (
     &    abs(s0h).gt.1.0d0/dble(myinum)
     &    ) then
        print*,''
        print*,"*** Warning in GENPHO: Minimum of horizontal beta-function not in origin ***"
        print*,''
      endif

      if (
     &    abs(s0v).gt.1.0d0/dble(myinum)
     &    ) then
        print*,''
        print*,"*** Warning in GENPHO: Minimum of vertical beta-function not in origin ***"
        print*,''
      endif

      if(
     &    abs(disp0).gt.1.0d-6
     &    .or.
     &    abs(ddisp0).gt.1.0d-6) then
        print*,''
        print*,"*** Warning in GENPHO: Disperion or it's derivative not zero ***"
        print*,''
      endif

      if (nsource.ne.1) then
        print*,''
        print*,"*** Warning in GENPHO: Number of sources not one ***"
        print*,''
      endif

      nedz=(nobsvz-mobsvz)/2
      nedy=(nobsvy-mobsvy)/2
      iobs=0
      jobs=0

      if (abs(genphophsh).eq.9999.0d0) then
        do ifrq=1,nfreq
          iobfr=icbrill+nobsv*(ifrq-1)
          rea(1:2)=(0.0d0,0.0d0)
          rea(3)=dcmplx(reaima(3,1,iobfr),reaima(3,2,iobfr))
          expsh=rea(3)/abs(rea(3))
          if (genphophsh.eq.-9999.0d0) expsh=expsh*cdexp(dcmplx(0.0d0,-pi1/2.0d0))
          do iy=1,nobsvy
            do iz=1,nobsvz
              iobs=iobs+1
              if(
     &          iz.le.nedz.or.iz.gt.nobsvz-nedz
     &          .or.
     &          iy.le.nedy.or.iy.gt.nobsvy-nedy
     &          ) cycle
              jobs=jobs+1
              obs(:,jobs)=obsv(:,iobs)
              iobfr=iobs+nobsv*(ifrq-1)
              jobfr=jobs+mobsv*(ifrq-1)
              rea(1:8)=dcmplx(reaima(1:8,1,iobfr),reaima(1:8,2,iobfr))/expsh
              reaima(1:8,1,iobfr)=dreal(rea)
              reaima(1:8,2,iobfr)=dimag(rea)
              arad(1:3,jobfr)=rea(1:3)
              arad(4:6,jobfr)=rea(6:8)
            enddo
          enddo
        enddo
      else
        expsh=cdexp(dcmplx(0.0d0,genphophsh))
        do ifrq=1,nfreq
          do iy=1,nobsvy
            do iz=1,nobsvz
              iobs=iobs+1
              if(
     &          iz.le.nedz.or.iz.gt.nobsvz-nedz
     &          .or.
     &          iy.le.nedy.or.iy.gt.nobsvy-nedy
     &          ) cycle
              jobs=jobs+1
              obs(:,jobs)=obsv(:,iobs)
              iobfr=iobs+nobsv*(ifrq-1)
              jobfr=jobs+mobsv*(ifrq-1)
              rea=dcmplx(reaima(1:8,1,iobfr),reaima(1:8,2,iobfr))/expsh
              reaima(1:8,1,iobfr)=dreal(rea)
              reaima(1:8,2,iobfr)=dimag(rea)
              arad(1:3,jobfr)=rea(1:3)
              arad(4:6,jobfr)=rea(6:8)
            enddo
          enddo
        enddo
      endif

      SPECNOR_SI= !merke/synchrotron_radiation.txt
     &  dmycur ! Strom
     &  /echarge1/hbar1*clight1/PI1*EPS01
     &  *banwid !BW

      sigz=sqrt(eps0h*beta0h)
      sigzp=sqrt(eps0h/beta0h)

      sigy=sqrt(eps0v*beta0v)
      sigyp=sqrt(eps0v/beta0v)

      ephomin=freq(1)
      ephomax=freq(nfreq)

      call urad_phase_amp_genpho(zi,yi,mobsvy,mobsvz,obs,
     &  moderan,nelecgenpho,noranone,npho,nfreq,
     &  ephomin,ephomax,
     &  sigz,sigzp,sigy,sigyp,
     &  ndimpho,photons,ndimele,electrons,
     &  arad,specnor_si)

      !allutil_break

      open(newunit=lunapho,file='ampgenpho.pho')
      open(newunit=lunaele,file='ampgenpho.elc')

      l=1
c      zmin=photons(2)
c      zmax=photons(ndimpho-ndimapho+2)
c      ymin=photons(3)
c      ymax=photons(ndimpho-ndimapho+3)

c      if (ipin.ne.0.and.(obsvz(1).eq.obsvz(nobsvz).or.obsvy(1).eq.obsvy(nobsvy))) then
c        kpin=0
c      else
c        kpin=ipin
c      endif

c      if (kpin.ne.0) then
c        zmin=obsvz(2)-obsvdz/10.
c        zmax=obsvz(nobsvz-1)+obsvdz/10.
c        ymin=obsvy(2)-obsvdy/10.
c        ymax=obsvy(nobsvy-1)+obsvdy/10.
c      endif

      do iel=1,nelecgenpho
        do i=1,ngam
          do iepho=1,nfreq
            z=photons(l+1)
            y=photons(l+2)
c            if(kpin.ne.0) then
c              if (z.lt.zmin.or.z.gt.zmax.or.y.lt.ymin.or.y.gt.ymax) then
c                l=l+ndimapho
c                cycle
c              endif
c            endif
            write(lunapho,*) i,iel,iepho,iefold,dmyenergy,g,
     &        photons(l:l),
     &        photons(l+1:l+4)*1000.0,
     &        photons(l+5:l+8)
            l=l+ndimapho
          enddo
        enddo
      enddo

      l=1
      do i=1,nelecgenpho
        write(lunaele,*) i,dmyenergy,g,
     &    electrons(l:l+1)*1000.0,
     &    electrons(l+2:l+3)*1000.0
        l=l+4
      enddo

      close(lunapho)
      close(lunaele)

      RETURN
      END
