*CMZ :  4.00/13 30/11/2021  11.08.07  by  Michael Scheer
*-- Author :    Michael Scheer   09/05/2012
*KEEP,waveg1.
      include 'waveg1.cmn'
*KEND.
      program wave_photons_main
      implicit none
*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      double precision :: x=0.0d0, y=0.0d0, z=0.0d0,
     &  veln(3)=[1.0d0,0.0d0,0.0d0],
     &  bx=0.0d0,by=1.30d0,bz=0.0d0, Ebeam=1.7,ds=0.001,
     &  dtim,dgamma,gamma

      integer :: i,isize=64,iseeds(64)

*KEEP,photon.
      include 'photon.cmn'
*KEND.

      gamma=Ebeam/emassg1
      nbing1=1000
      eecmaxg1=5.0
      iphmode=-1
      dtim=ds/clight1

      call util_random_init(isize,iseeds)

      open(unit=99,file='photons_wave.dat')

      qfmean=qfmean/nqfphotons*1.0e6 !KeV
      qfrms=sqrt((qfrms/nqfphotons-qfmean**2))*1.0e6 !KeV

      do i=1,100000
        call photon_simple(x,y,z,veln,gamma,bx,by,bz,dgamma,dtim)
        if (norm2(dpphoton).gt.0.0d0) then
c          print '(I10,3e15.6)',i,dpphoton*1.0e6 !KeV
          write(99,'(5e15.6)') dpphoton*1.0e6,qfmean,qfrms !KeV
        endif
      enddo

      close(99)

      call util_random_get_seed(isize,iseeds)

      open(unit=99,file='util.seeds')

      write(99,*)isize,"-9999 ! written by util_random_main.f"
      do i=1,isize
        write(99,*)i,iseeds(i)
      enddo
      close(99)

      end

! +PATCH,//WAVE/MAIN
! +DECK,wave_photons_main.

      subroutine photon_simple(x,y,z,veln,gamma,bx,by,bz,dgamma,dtim)

      use wave_g1

      implicit none

      integer ical,i,ncy,ieof,iwarn

      double precision x,y,z,veln(3),bmag(3),bx,by,bz,ebeam,elmom,gamma,
     &  bparn,bper(3),bpern,epho,eec,bpervn(3),
     &  dgamma,b2per,
     &  pdum,dtim,ec,photons,de,deecg1,eecg1,g1,yrnint10,sigv,dum

      real rnrn(2),hrndm1m,xran(1)
      real*8 fill(100)
      double precision, dimension (:), allocatable ::
     &  xrn,yrn,yrnint,coef,work1,work2,work3,work4,
     &  coefpsi,eecpsilog,sigpsi,cylog

*KEEP,photon.
      include 'photon.cmn'
*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      data ical/0/,iwarn/0/,ncy/0/

      save

      sigv=0.0d0
      nbing1=max(nbing1,10)

      if (ical.eq.0) then

        ncy=11 !27.4.2020

        allocate(coefpsi(max(nbing1,ncy)))
        allocate(eecpsilog(max(nbing1,ncy)))
        allocate(cylog(max(nbing1,ncy)))
        allocate(sigpsi(max(nbing1,ncy)))

        do i=1,ncy
          eecpsilog(i)=eecpsilog1(i)
          sigpsi(i)=sigpsi1(i)
          cylog(i)=cylog1(i)
        enddo

        eecpsilog(1:ncy)=log(eecpsilog(1:ncy))
        cylog(1:ncy)=log(cylog(1:ncy)/1000.0d0) !mrad -> rad

        qfrms=0.0d0
        qfmean=0.0d0
        nqfphotons=0

        if (nbing1.le.2) then
          write(6,*)'*** Error in photon: NBING1 must be greater then 1'
          stop
        endif

        allocate(xrn(max(nbing1,ncy)))
        allocate(yrn(max(nbing1,ncy)))
        allocate(yrnint(max(nbing1,ncy)))
        allocate(coef(max(nbing1,ncy)))
        allocate(work1(max(nbing1,ncy)))
        allocate(work2(max(nbing1,ncy)))
        allocate(work3(max(nbing1,ncy)))
        allocate(work4(max(nbing1,ncy)))

        deecg1=eecmaxg1/(nbing1-1)

        eecg1=0.0d0
        do i=1,10
          eecg1=eecg1+deecg1/10.0d0
          call util_g1_static(eecg1,g1)
          xrn(i)=eecg1
          yrn(i)=g1/eecg1
        enddo

        yrnint(1)=0.0d0
        do i=2,10
          yrnint(i)=yrnint(i-1)+(yrn(i)+yrn(i-1))/2.0d0*(xrn(i)-xrn(i-1))
        enddo
        yrnint10=yrnint(10)

        eecg1=0.0d0
        do i=10,nbing1
          eecg1=eecg1+deecg1
          call util_g1_static(eecg1,g1)
          xrn(i)=eecg1
          yrn(i)=g1/eecg1
        enddo

        call util_spline_running_integral(
     &    xrn(10:nbing1),yrn(10:nbing1),nbing1-10+1,yrnint(10:nbing1),
     &    coef,work1,work2,work3,work4)

        yrnint(10)=yrnint10
        yrnint(11:nbing1)=yrnint(11:nbing1)+yrnint10
        yrnint=yrnint/yrnint(nbing1)

        do i=2,nbing1
          if (yrnint(i).le.yrnint(i-1)) then
            stop '*** Error in photon: Bad integration of G1 ***'
          endif
        enddo

        call util_spline_coef(
     &    yrnint,xrn,nbing1,0.0d0,0.0d0,
     &    coef,work1,work2,work3,work4)

        if (iphmode.lt.0) then
          call util_spline_coef(
     &      eecpsilog,cylog,ncy,0.0d0,0.0d0,
     &      coefpsi,work1,work2,work3,work4)
        endif

        !dgamma=-pdum*gamma**2*b**2*dt
        pdum=cgam1/2.0d0/pi1*clight1*(clight1/1.0d9)**2*emassg1

        ical=1

      endif !ical

      bmag(1)=bx
      bmag(2)=by
      bmag(3)=bz

      elmom=emassg1*dsqrt((gamma-1.0d0)*(gamma+1.0d0)) !GeV
      ebeam=emassg1*gamma !GeV

      bparn=(bmag(1)*veln(1)+bmag(2)*veln(2)+bmag(3)*veln(3))
      bper=bmag-bparn*veln
      b2per=bper(1)**2+bper(2)**2+bper(3)**2
      bpern=sqrt(b2per)
      bpervn=bper/bpern

      ec=ecdipkev1*bpern*ebeam**2*1.0d-6 !GeV

      !dgamma = pdum * gamma**2 * b2per * dtim
      !dN = 15*sqrt(3)/8 * dE/Ec = 3.2476 * de/ec

      de=pdum*b2per*gamma*ebeam*dtim !GeV

      if (ec.ne.0.0d0) then
        photons=3.2476d0*de/ec !number of photons
      else
        photons=0.0d0
      endif

      if (photons.ge.1.0d0.and.iwarn.eq.0) then
        write(6,*)'*** Warning in PHOTON: Step size to large, ***'
        write(6,*)'*** i.e. probabilty to generate photon is greater than one! ***'
        iwarn=1
      endif

      call util_random(2,rnrn)

      if(rnrn(1).le.photons) then

        call util_spline_inter(yrnint,xrn,coef,nbing1,
     &    dble(rnrn(2)),eec,-1)
        if (eec.lt.0.0d0) then
          print*,'*** Warning in PHOTON: Negative photon energy occured ***'
          print*,'rnrn:',rnrn
          print*,'setting Epho/Ec = 1.e-6'
          eec=1.0d-6
        endif
        if (iphmode.lt.0) then
          call util_spline_inter(eecpsilog,cylog,coefpsi,ncy,log(eec),
     &      sigv,-1)
          sigv=0.408*exp(sigv)/ebeam !rad
          call util_random_gauss(1,xran)
          sigv=xran(1)*sigv
        else
          sigv=0.0d0
        endif !iphmode

        epho=eec*ec
        dgamma=-epho/ebeam*gamma
        dpphoton=-(veln+sigv*bpervn)*epho

        qfmean=qfmean+epho
        qfrms=qfrms+epho**2

        nqfphotons=nqfphotons+1

      else
        dpphoton=0.0d0
        dgamma=0.0d0
      endif !(rnrn.le.wrad)

      return
      end

      include 'util_random_get_seed.f'
      include 'util_random_init.f'
      include 'util_random.f'
