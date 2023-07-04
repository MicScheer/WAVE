*CMZ :  4.00/13 16/11/2021  21.29.08  by  Michael Scheer
*CMZ :  3.05/05 09/07/2018  13.41.08  by  Michael Scheer
*CMZ :  3.05/04 05/07/2018  12.40.34  by  Michael Scheer
*CMZ :  3.03/04 01/12/2017  14.34.43  by  Michael Scheer
*CMZ :  3.02/04 13/03/2015  10.30.49  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  15.40.59  by  Michael Scheer
*CMZ :  2.68/05 07/09/2012  11.23.42  by  Michael Scheer
*CMZ :  2.68/04 03/09/2012  11.49.18  by  Michael Scheer
*CMZ :  2.68/03 31/08/2012  09.22.41  by  Michael Scheer
*-- Author : Michael Scheer
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
      program urad_main

      implicit none

      integer nstepp,nphp
      parameter (nstepp=100000,nphp=10001)

      double precision
     &  gammai,dgamtot,powden,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xf,yf,zf,efx,efy,efz,
     &  xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,ds,
     &  traxyz(14,nstepp),
     &  xobsv,yobsv,zobsv,
     &  phelow,phehig,phener(nphp),stokes(4,nphp)

      double precision energy,omegac

      complex*16 aradx(nphp),arady(nphp),aradz(nphp)

      real xpl(nstepp),ypl(nstepp),zpl(nstepp)

      integer ::
     &  nstep,nthstep,ndim,ifreq,istep,icharge=-1,
     &  nphener,ieneloss,ivelofield,istatus

      character(2048) cline

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEEP,URADNAM.
      namelist/uradn/
     &  energy,xobsv,yobsv,zobsv,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xf,yf,zf,efx,efy,efz,ds,
     &  phelow,phehig,
     &  nstep,nphener,ieneloss,ivelofield,nthstep
*KEEP,URADCOM.
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
*KEEP,random.
      include 'random.cmn'
*KEND.

      save
     &  traxyz,powden,stokes,aradx,arady,aradz,
     &  xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,dgamtot,phener,
     &  istatus

      ndim=nstepp

      open(unit=99,file='urad.nam',status='old')
      read(99,uradn)
      close(99)

      kbmapu=kbmap
      kmagsequ=kmagseq

      irnmode=0
      irnsize=64
      call util_random_init(irnsize,irnseed)

      kbmapu=kbmap
      kmagsequ=kmagseq
      gammai=energy/emassg1

      omegac=1.5d0*gamma1**3*clight1/rho1

      call zeit(6)

      call urad(icharge
     &  ,gammai,dgamtot,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xf,yf,zf,efx,efy,efz,
     &  xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,ds
     &  ,nthstep,nstep,ndim
     &  ,traxyz
     &  ,xobsv,yobsv,zobsv
     &  ,phelow,phehig
     &  ,nphener
     &  ,phener
     &  ,aradx,arady,aradz,
     &  stokes,powden,ieneloss,ivelofield
     &  ,istatus
     &  )

      if (istatus.ne.0) then
        print*,"*** Error in utest: Bad return from utest, istatus:",istatus
        print*,"istatus = -2: Dimension of traxyz exceeded!"
      endif

      call zeit(6)

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

      open(unit=99,file='urad.nam',status='old')
      do while (.true.)
        read(99,*,end=99)cline
        print*,trim(cline)
      enddo
99    close(99)

      print*,' '
      print*,'xexit, yexit, zexit:',sngl(xexit),sngl(yexit),sngl(zexit)
      print*,' '
      print*,'vnxex, vnyex, vnzex:',sngl(vnxex),sngl(vnyex),sngl(vnzex)
      print*,' '
      print*,'rel. energy-loss:',dgamtot/gammai
      print*,'power density:',powden

      print*,'istatus:',istatus

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
      endif

      call urad_plot_traj(14,nstepp,nstep,traxyz,xpl,ypl,zpl)
      call urad_plot_b(14,nstepp,nstep,traxyz,xpl,ypl,zpl)
      call urad_plot_stokes(nphener,phener,stokes,xpl,ypl)

      end
