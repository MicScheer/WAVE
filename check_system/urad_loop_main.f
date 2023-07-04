*CMZ :  3.05/04 28/06/2018  09.44.26  by  Michael Scheer
*CMZ :  3.04/00 24/01/2018  13.25.20  by  Michael Scheer
*CMZ :  3.03/04 18/12/2017  16.17.42  by  Michael Scheer
*-- Author :    Michael Scheer   11/10/2017

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

      program urad_loop_main

      implicit none

      integer nstepp,nphp
      parameter (nstepp=10000000,nphp=10001)

      double precision
     &  gammai,dgamtot,powden,beta,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xelec0,yelec0,zelec0,
     &  vxelec0,vyelec0,vzelec0,
     &  xf,yf,zf,efx,efy,efz,
     &  vxexit,vyexit,vzexit,
     &  xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,ds,
     &  traxyz(14,nstepp),
     &  xobsv,yobsv,zobsv,
     &  phelow,phehig,phener(nphp),stokes(4,nphp)
     &  ,zi,zpi,yi,ypi,zpf,ypf
     &  ,ze,zpe,ye,ype,xib,yib,zib,vxelecb,vyelecb,vzelecb,
     &  xcoi,ycoi,zcoi,vxcoi,vycoi,vzcoi,
     &  xcof,ycof,zcof,vxcof,vycof,vzcof,
     &  vxi,vyi,vzi,vxf,vyf,vzf,xi,
     &  vnxcoi,vnycoi,vnzcoi,vnxcof,vnycof,vnzcof,
     &  emiti,emitf,det,dmati(4,4),dmatf(4,4),zcogf,ycogf,zpcogf,ypcogf,
     &  avezi,aveyi,avezpi,aveypi,avezf,aveyf,avezpf,aveypf,sexit

      double precision :: scalmat=1.0d6
      double precision, dimension (:,:), allocatable :: zzpyypi,zzpyypf

      double precision
     &  energy,vn

      complex*16
     &  aradx(nphp),arady(nphp),aradz(nphp)

      integer iutil_fexist,
     &  nstep,ndim,ifreq,istep,ifail,
     &  nphener,ieneloss,ivelofield,istatus,isample,
     &  moop,nloop,loop,lulo,icorrco,icharge,igenfun,itrackback

      namelist/uradn/
     &  energy,xobsv,yobsv,zobsv,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xf,yf,zf,efx,efy,efz,ds,
     &  phelow,phehig,
     &  nstep,nphener,ieneloss,ivelofield,icharge,igenfun,itrackback,
     &  isample

*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,uradcom.
      integer nfour,nfourwls,ifour0,iprntf,maxfoumagp,nfoumags,nfoumagcp
      parameter (maxfoumagp=100,nfoumagcp=2**10)

      double precision
     &  xfoubounds(5,maxfoumagp),foumags(nfoumagcp/2+3,maxfoumagp)
     &  ,fouentr,fouexit,xshbfour

c      character(2048) chfoumags(maxfoumagp)

      integer kmonopole,intpolbmap,kbmap,kmagseq,kbmapu,kmagsequ

      double precision xlenfour,xbhomf,dbhomf,emom
      character(2048) fmagseq

      common/cfourier/
     &  xlenfour,xbhomf,dbhomf,emom
     &  ,xfoubounds,foumags
     &  ,fouentr,fouexit,xshbfour
     &  ,kmagseq,fmagseq
     &  ,nfour,nfourwls,ifour0,iprntf,nfoumags,kmonopole

      integer iahwfour,nhhalbasy

      double precision b0scale,hshift,vshift

      integer khalbasy

      common/bscalec/ b0scale,hshift,vshift,intpolbmap

      double precision b0halbasy,
     &  xlhalbasy,ylhalbasy,zlhalbasy,
     &  xkhalbasy,ykhalbasy,zkhalbasy,
     &  fasym,ahwpol,rhalbasy,xcenhal,hhalbasy,pkhalbasy,
     &  ugamma,uenergy,ucur

      common/uradhalbasym/
     &  b0halbasy,
     &  xlhalbasy,ylhalbasy,zlhalbasy,
     &  xkhalbasy,ykhalbasy,zkhalbasy,
     &  fasym,ahwpol,rhalbasy,xcenhal,hhalbasy,pkhalbasy,
     &  ugamma,uenergy,ucur,
     &  khalbasy,iahwfour

      integer nmgsqp,mmag,irfilf,iwfilf
      parameter(nmgsqp=10000)

      character(3) ctyp(nmgsqp)

      double precision pmag(13,nmgsqp),coz(2,nmgsqp+1)
     &  ,corr(nmgsqp),uebounds(2,nmgsqp),
     &  dibounds(2,nmgsqp),
     &  dhbounds(2,nmgsqp),
     &  qfbounds(2,nmgsqp),
     &  sxbounds(2,nmgsqp),
     &  bmsqbounds(2)

      common/uradmgsqc/pmag,coz,corr,dibounds,dhbounds,
     &  uebounds,qfbounds,sxbounds,bmsqbounds,mmag,irfilf,iwfilf,
     &  kbmapu,kmagsequ,
     &  ctyp

      namelist/ufield/nfour,nfourwls,ifour0,xlenfour,dbhomf,iprntf,irfilf,
     &  b0halbasy,xlhalbasy,ylhalbasy,hhalbasy,pkhalbasy,
     &  zlhalbasy,fasym,ahwpol,iahwfour,xcenhal,nhhalbasy,
     &  b0scale,hshift,vshift,kmonopole,intpolbmap,kbmap,kmagseq,fmagseq
*KEND.

      save

*KEEP,phycon1.
      include 'phycon1.cmn'
*KEND.

      ndim=nstepp
      emom=emasse1*dsqrt((gammai-1.0d0)*(gammai+1.0d0)) !eV

      open(unit=99,file='urad.nam',status='old')
      read(99,uradn)
      read(99,uradfoun)
      close(99)

      if (isample.eq.0) isample=1

      if (efx.ne.1.0d0) then
        print*
        print*,"*** Warning in urad_loop_main.f: efx .ne. 0 not implemented ***"
        print*
      endif

      gammai=energy/emassg1
      beta=dsqrt((1.0d0-1.0d0/gammai)*(1.0d0+1.0d0/gammai))

      xcoi=xelec
      ycoi=yelec
      zcoi=zelec

      vn=sqrt(vxelec**2+vyelec**2+vzelec**2)
      vnxcoi=vxelec/vn
      vnycoi=vyelec/vn
      vnzcoi=vzelec/vn
      vxcoi=vxelec/vn*beta*clight1
      vycoi=vyelec/vn*beta*clight1
      vzcoi=vzelec/vn*beta*clight1

      if (icharge.gt.0) icharge=1
      if (icharge.le.0) icharge=-1

      call zeit(6)

      print*
      print*,"First call urad to get reference orbit"
c      print*,"First call to utrack to get reference orbit"
      print*

c      call utrack(icharge,
c     &  gammai,dgamtot,
c     &  xcoi,ycoi,zcoi,vxcoi,vycoi,vzcoi,
c     &  xf,yf,zf,efx,efy,efz,
c     &  xcof,ycof,zcof,vnxex,vnyex,vnzex,texit,ds,
c     &  ieneloss,ivelofield
c     &  ,istatus)

      call urad(icharge,
     &  gammai,dgamtot,
     &  xcoi,ycoi,zcoi,vxcoi,vycoi,vzcoi,
     &  xf,yf,zf,efx,efy,efz,
     &  xcof,ycof,zcof,vnxex,vnyex,vnzex,texit,ds,
     &  nstep,ndim,traxyz,
     &  xobsv,yobsv,zobsv,
     &  phelow,phehig,nphener,
     &  phener,aradx,arady,aradz,stokes,powden,ieneloss,ivelofield
     &  ,istatus)

      vnxcof=vnxex
      vnycof=vnyex
      vnzcof=vnzex

      vxcof=vnxcof*beta*clight1
      vycof=vnycof*beta*clight1
      vzcof=vnzcof*beta*clight1

      zcogf=zcof
      ycogf=ycof
      zpcogf=vzcof/vxcof
      ypcogf=vycof/vxcof

      if (igenfun.ne.0) then

        call urad_idtrmshgf(igenfun,nfour,hshift,vshift,
     &    xcoi,ycoi,zcoi,vnxcoi,vnycoi,vnzcoi,
     &    xcof,ycof,zcof,vnxcof,vnycof,vnzcof,
     &    icharge,gammai,
     &    xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &    xexit,yexit,zexit,vxexit,vyexit,vzexit,sexit,
     &    istatus)

        zcogf=zexit
        ycogf=yexit
        zpcogf=vzexit/vxexit
        ypcogf=vyexit/vxexit

      endif

      open(unit=99,file='urad_reference_orbit.dat',recl=256)
        write(99,*)xcoi,ycoi,zcoi,vxcoi,vycoi,vzcoi
        write(99,*)xcof,ycof,zcof,vxcof,vycof,vzcof
      close(99)

      open(unit=99,file='urad_traxyz.dat',recl=256)
      do istep=1,nstep
        write(99,*)sngl(traxyz(1:14,istep))
      enddo
      close(99)

      open(unit=99,file='urad_stokes.dat',recl=256)
      do ifreq=1,nphener
        write(99,*)sngl(phener(ifreq)),sngl(stokes(1:4,ifreq))
      enddo
      close(99)

      open(unit=99,file='urad_amplitude.dat',recl=256)
      do ifreq=1,nphener
        write(99,*)sngl(phener(ifreq))
     &    ,sngl(dreal(aradx(ifreq))),sngl(dimag(aradx(ifreq)))
     &    ,sngl(dreal(arady(ifreq))),sngl(dimag(arady(ifreq)))
     &    ,sngl(dreal(aradz(ifreq))),sngl(dimag(aradz(ifreq)))
      enddo
      close(99)

      open(newunit=lulo,file='loop.out')
      write(lulo,'(a)')"* loop,xi,yi,zi,vxi,vyi,vzi,xf,yf,zf,vxf,vyf,vzf,istat"

      call zeit(6)

      open(unit=99,file='loop.in',status='old')
      call util_skip_comment(99)
      read(99,*)nloop,zi,zpi,yi,ypi,icorrco
      close(99)

      xelec=xcoi
      yelec=yi
      zelec=zi
      vxelec=beta*clight1/sqrt(1.0d0+(ypi**2+zpi**2))
      vyelec=vxelec*ypi
      vzelec=vxelec*zpi

      print*
      print*,"Starting tracking loop"
      print*

      allocate(zzpyypi(4,nloop))
      allocate(zzpyypf(4,nloop))

      moop=0
      do loop=1,nloop

        xi=xelec
        yi=yelec
        zi=zelec

        ypi=vyelec/vxelec
        zpi=vzelec/vxelec

        vxi=beta*clight1/sqrt(1.0d0+(ypi**2+zpi**2))
        vyi=vxi*ypi
        vzi=vxi*zpi

        if (igenfun.eq.0) then

          call utrack(icharge,
     &      gammai,dgamtot,
     &      xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &      xf,yf,zf,efx,efy,efz,
     &      xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,ds,
     &      ieneloss,ivelofield
     &      ,istatus)

c          call urad(icharge,
c     &      gammai,dgamtot,
c     &      xelec,yelec,zelec,vxelec,vyelec,vzelec,
c     &      xf,yf,zf,efx,efy,efz,
c     &      xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,ds,
c     &      nstep,ndim,traxyz,
c     &      xobsv,yobsv,zobsv,
c     &      phelow,phehig,nphener,
c     &      phener,aradx,arady,aradz,stokes,powden,ieneloss,ivelofield
c     &      ,istatus)

          vxf=vnxex*beta*clight1
          vyf=vnyex*beta*clight1
          vzf=vnzex*beta*clight1

          if (icorrco.ne.0) then
            yelec=yexit-ycof
            zelec=zexit-zcof
            zpf=vzf/vxf-vzcof/vxcof
            ypf=vyf/vxf-vycof/vxcof
            vxelec=vxf
            vyelec=vxf*ypf
            vzelec=vxf*zpf
          else
            yelec=yexit
            zelec=zexit
            vxelec=vxf
            vyelec=vyf
            vzelec=vzf
          endif

        else

          vn=sqrt(vxelec**2+(vyelec**2+vzelec**2))

          vxelec=vxelec/vn*clight1*beta
          vyelec=vyelec/vn*clight1*beta
          vzelec=vzelec/vn*clight1*beta

          call urad_idtrmshgf(igenfun,nfour,hshift,vshift,
     &      xcoi,ycoi,zcoi,vnxcoi,vnycoi,vnzcoi,
     &      xcof,ycof,zcof,vnxcof,vnycof,vnzcof,
     &      icharge,gammai,
     &      xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &      xexit,yexit,zexit,vxexit,vyexit,vzexit,sexit,
     &      istatus)

          vn=sqrt(vxexit**2+(vyexit**2+vzexit**2))

          vnxex=vxexit/vn
          vnyex=vyexit/vn
          vnzex=vzexit/vn

          vxf=vnxex*beta*clight1
          vyf=vnyex*beta*clight1
          vzf=vnzex*beta*clight1

          if (icorrco.ne.0) then
            yelec=yexit-ycogf
            zelec=zexit-zcogf
            zpf=vzf/vxf-zpcogf
            ypf=vyf/vxf-ypcogf
            vxelec=vxf
            vyelec=vxf*ypf
            vzelec=vxf*zpf
          else
            yelec=yexit
            zelec=zexit
            vxelec=vnxex
            vyelec=vnyex
            vzelec=vnzex
          endif

        endif

        if (mod(loop,isample).eq.0) then
          write(lulo,*)loop,
     &      xi,yi,zi,vxi,vyi,vzi,
     &      xexit,yexit,zexit,vxf,vyf,vzf,istatus
          if (iutil_fexist("urad.cont").eq.0) then
            print*,"*** Exiting loop since urad.cont is missing ***"
            exit
          endif
        endif

        if (istatus.ne.0) then
          print*,"*** Bad return from urad_idtrmshgf, exiting loop ***"
          exit
        endif

        moop=moop+1

        zzpyypi(1,moop)=zi-zcoi
        zzpyypi(2,moop)=vzi/vxi-vzcoi/vxcoi
        zzpyypi(3,moop)=yi-ycoi
        zzpyypi(4,moop)=vyi/vxi-vycoi/vxcoi

        zzpyypf(1,moop)=zexit-zcogf
        zzpyypf(2,moop)=vzf/vxf-zpcogf
        zzpyypf(3,moop)=yexit-ycogf
        zzpyypf(4,moop)=vyf/vxf-ypcogf

      enddo

      call zeit(6)

      print*
      print*,"Tracking loop done"
      print*

      avezi=0.0d0
      avezpi=0.0d0
      aveyi=0.0d0
      aveypi=0.0d0
      avezf=0.0d0
      avezpf=0.0d0
      aveyf=0.0d0
      aveypf=0.0d0
      dmati=0.0d0
      dmatf=0.0d0

      do loop=1,moop
        avezi=avezi+zzpyypi(1,loop)
        avezpi=avezpi+zzpyypi(2,loop)
        aveyi=aveyi+zzpyypi(3,loop)
        aveypi=aveypi+zzpyypi(4,loop)
        avezf=avezf+zzpyypf(1,loop)
        avezpf=avezpf+zzpyypf(2,loop)
        aveyf=aveyf+zzpyypf(3,loop)
        aveypf=aveypf+zzpyypf(4,loop)
      enddo

      avezi=avezi/moop
      avezpi=avezpi/moop
      aveyi=aveyi/moop
      aveypi=aveypi/moop

      avezf=avezf/moop
      avezpf=avezpf/moop
      aveyf=aveyf/moop
      aveypf=aveypf/moop

      do loop=1,moop
        dmati(1,1)=dmati(1,1)+(zzpyypi(1,loop)-avezi)*(zzpyypi(1,loop)-avezi)
        dmati(1,2)=dmati(1,2)+(zzpyypi(1,loop)-avezi)*(zzpyypi(2,loop)-avezpi)
        dmati(1,3)=dmati(1,3)+(zzpyypi(1,loop)-avezi)*(zzpyypi(3,loop)-aveyi)
        dmati(1,4)=dmati(1,4)+(zzpyypi(1,loop)-avezi)*(zzpyypi(4,loop)-aveypi)
        dmati(2,1)=dmati(1,2)
        dmati(2,2)=dmati(2,2)+(zzpyypi(2,loop)-avezpi)*(zzpyypi(2,loop)-avezpi)
        dmati(2,3)=dmati(2,3)+(zzpyypi(2,loop)-avezpi)*(zzpyypi(3,loop)-aveyi)
        dmati(2,4)=dmati(2,4)+(zzpyypi(2,loop)-avezpi)*(zzpyypi(4,loop)-aveypi)
        dmati(3,1)=dmati(1,3)
        dmati(3,2)=dmati(2,3)
        dmati(3,3)=dmati(3,3)+(zzpyypi(3,loop)-aveyi)*(zzpyypi(3,loop)-aveyi)
        dmati(3,4)=dmati(3,4)+(zzpyypi(3,loop)-aveyi)*(zzpyypi(4,loop)-aveypi)
        dmati(4,1)=dmati(1,4)
        dmati(4,2)=dmati(2,4)
        dmati(4,3)=dmati(3,4)
        dmati(4,4)=dmati(4,4)+(zzpyypi(4,loop)-aveypi)*(zzpyypi(4,loop)-aveypi)

        dmatf(1,1)=dmatf(1,1)+(zzpyypf(1,loop)-avezf)*(zzpyypf(1,loop)-avezf)
        dmatf(1,2)=dmatf(1,2)+(zzpyypf(1,loop)-avezf)*(zzpyypf(2,loop)-avezpf)
        dmatf(1,3)=dmatf(1,3)+(zzpyypf(1,loop)-avezf)*(zzpyypf(3,loop)-aveyf)
        dmatf(1,4)=dmatf(1,4)+(zzpyypf(1,loop)-avezf)*(zzpyypf(4,loop)-aveypf)
        dmatf(2,1)=dmatf(1,2)
        dmatf(2,2)=dmatf(2,2)+(zzpyypf(2,loop)-avezpf)*(zzpyypf(2,loop)-avezpf)
        dmatf(2,3)=dmatf(2,3)+(zzpyypf(2,loop)-avezpf)*(zzpyypf(3,loop)-aveyf)
        dmatf(2,4)=dmatf(2,4)+(zzpyypf(2,loop)-avezpf)*(zzpyypf(4,loop)-aveypf)
        dmatf(3,1)=dmatf(1,3)
        dmatf(3,2)=dmatf(2,3)
        dmatf(3,3)=dmatf(3,3)+(zzpyypf(3,loop)-aveyf)*(zzpyypf(3,loop)-aveyf)
        dmatf(3,4)=dmatf(3,4)+(zzpyypf(3,loop)-aveyf)*(zzpyypf(4,loop)-aveypf)
        dmatf(4,1)=dmatf(1,4)
        dmatf(4,2)=dmatf(2,4)
        dmatf(4,3)=dmatf(3,4)
        dmatf(4,4)=dmatf(4,4)+(zzpyypf(4,loop)-aveypf)*(zzpyypf(4,loop)-aveypf)
      enddo

      dmati=dmati/moop*scalmat
      dmatf=dmatf/moop*scalmat

      call util_determinante(4,dmati,det,ifail)
      if (ifail.ne.0.or.det.lt.0) then
        print*,"Deteminante negative of bad return from util_determinante:",det,ifail
        emiti=-9999.
      else
        emiti=16.0d0*sqrt(det/scalmat**4)
      endif

      call util_determinante(4,dmatf,det,ifail)
      if (ifail.ne.0.or.det.lt.0) then
        print*,"Deteminante negative of bad return from util_determinante:",det,ifail
        emitf=-9999.
      else
        emitf=16.0d0*sqrt(det/scalmat**4)
      endif

c back tracking

      if (igenfun.eq.0.and.itrackback.ne.0) then
        stop "*** Backtracking disabled ***"
        print*
        print*,"Starting back tracking loop"
        print*
        call zeit(6)
        print*
        print*,"Back tracking loop done"
        print*
      endif !igenfun

      close(lulo)

c      call system('cat urad.nam')
      print*,' '
      print*,'xexit, yexit, zexit:',sngl(xexit),sngl(yexit),sngl(zexit)
      print*,' '
      print*,'vnxex, vnyex, vnzex:',sngl(vnxex),sngl(vnyex),sngl(vnzex)
      print*,' '
      print*,'rel. energy-loss:',dgamtot/gammai
      print*,'power density:',powden

      print*
      print*,"Sigma matrix at entrance:"
      print*
      print '(4e15.5)',dmati(1,1:4)
      print '(4e15.5)',dmati(2,1:4)
      print '(4e15.5)',dmati(3,1:4)
      print '(4e15.5)',dmati(4,1:4)
      print*

      print*
      print*,"Sigma matrix at exit:"
      print*
      print '(4e15.5)',dmatf(1,1:4)
      print '(4e15.5)',dmatf(2,1:4)
      print '(4e15.5)',dmatf(3,1:4)
      print '(4e15.5)',dmatf(4,1:4)
      print*
      print*
      print*,"emittance at entrance and exit, 1-ratio:",
     &  sngl(emiti),sngl(emitf),
     &  sngl(1.0d0-emiti/emitf)

      print*
      print*,'istatus:',istatus
      print*

      if (istatus.eq.-1) then
        print*,'*** Error in urad: Zero gamma or velocity vector ***'
      else if (istatus.eq.-2) then
        print*,'***  Error in urad: Dimension of traxyz exceeded ***'
      else if (istatus.eq.-3) then
        print*,'***  Error in urad: Bad value of ivelofield, must be 0, 1, or 2 ***'
      else if (istatus.eq.999) then
        print*,'***  Error in uradbmap:  Too few data ***'
      else if (istatus.eq.1001) then
        print*,'***  Error in uradbmap:  Out of range in x ***'
      else if (istatus.eq.1002) then
        print*,'***  Error in uradbmap:  Out of range in y ***'
      else if (istatus.eq.1003) then
        print*,'***  Error in uradbmap:  Out of range in z ***'
      else if (istatus.eq.-11) then
        print*,'***  Instable particle in urad_idtrmshgf ***'
      endif

      stop
      end

      include 'urad.f'
      include 'uradfield.f'
      include 'uradrndm.f'
      include 'uradbmap.f'
      include 'urad_idtrmshgf.f'
      include 'uradbmagseq.f'
      include 'uradbmagseqc.f'
      include 'uradbdi.f'
      include 'uradbdh.f'
c      include 'utrack.f'
