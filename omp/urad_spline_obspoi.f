*CMZ :          30/12/2023  16.11.21  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  09.53.09  by  Michael Scheer
*-- Author :    Michael Scheer   02/12/2023
      subroutine urad_spline_obspoi(
     &  icharge,current,park,nper,perl,
     &  gammai,dgamtot,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xf,yf,zf,efxn,efyn,efzn,
     &  xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,dsi,
     &  nlpoi,traxyz,phase0,
     &  xobsv,yobsv,zobsv,phelow,phehig,
     &  nphener,phener,aradex,aradey,aradez,aradbx,aradby,aradbz,stokes,powden,
     &  ieneloss,ivelofield,
     &  istatus,ith,banwid,modewave,iobsv
     &  )
      implicit none

      double complex :: ay,az,zi=(0.0d0,1.0d0),expom,
     &  aradex(nphener),aradey(nphener),aradez(nphener),
     &  aradbx(nphener),aradby(nphener),aradbz(nphener),
     &  apolh,apolr,apoll,apol45,
     &  ampex,ampey,ampez,ampbx,ampby,ampbz,
     &  zone=(1.0d0,0.0d0)

      double precision
     &  current,gammai,dgamtot,charge,
     &  xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &  xf,yf,zf,efxn,efyn,efzn,
     &  xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,ds,dsi,
     &  phase0,xobsv,yobsv,zobsv,phelow,phehig,
     &  phener(nphener),stok1,stok2,stok3,stok4,
     &  stokes(4,nphener),powden,
     &  banwid,pownor,rspn,specnor,yp1,ypn,
     &  b3,bet1n,bpx,bpy,bpz,br2,br4,dom1,dom2,dphase,dum11,dum3,gamma,perl

      integer :: modewave,kfreq,nlpoi,nallo=0,iobsv,i,
     &  icharge,nthstep,nstep,ieneloss,ivelofield,istatus,ith,nphener

      double precision :: park,xlen,beta,betx,bety,betz,om,wlen1,eharm1,
     &  tharm1,txlen,Tu,tdev,zp,t,dt,dt2,tk,x,y,z,s,xk,xkx,dx,v0,yp,dom,
     &  rnnb(3),rnx,rny,rnz,r,azr,azi,ayr,ayi,eyr,eyi,ezr,ezi,h2,dum,phase,
     &  px,py,pz,r1,rnbx,rnby,rnbz,rn2,rn4,rx,ry,rz,vx2,vy2,vz2,vxp,vyp,vzp,
     &  x2,y2,z2,rarg5(5),rnr2,rnr4,vx12,vy12,vz12,acc,bx,by,bz,
     &  ct,st,betz0,zamp

      double precision traxyz(14,nlpoi)
      double complex, dimension (:,:), allocatable :: expo
      double precision, dimension (:,:), allocatable :: amp
      double precision, dimension (:,:), allocatable :: rarg
      double precision, dimension (:), allocatable ::
     &  tspl,ayrspl,ayispl,azrspl,azispl,w1,w2,w3,w4,w5,y2p

      integer :: nper,nx,ix,kx,luno,luna,ifr,ifail,k,idebug=0,ical=0

      save

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,observ.
      include 'observ.cmn'
*KEND.
c+seq,freqs.
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      ical=ical+1

      if (idebug.ne.0) then
c        print*,"Entered urad_spline_obspoi, ical:",ical,idebug
      endif

      if (idebug.gt.1.and.yobsv.eq.0.0d0.and.zobsv.eq.0.0d0) then
        open(newunit=luno,file='urad_spline.out',recl=1024)
        open(newunit=luna,file='urad_spline.amp',recl=1024)
      endif

      aradex=dcmplx((0.0d0,0.0d0))

      gamma=gammai
      beta=dsqrt((1.0d0-1.0d0/gammai)*(1.0d0+1.0d0/gammai))
      v0=beta*clight1

      xlen=xf-xelec
      wlen1=(1.0d0+park**2/2.0d0)/2.0d0/gammai**2*xlen
      eharm1=wtoe1/(wlen1*1.0d9)
      tharm1=wlen1/clight1

      tdev=(xlen+wlen1)/clight1
      tu=tdev/dble(nper)

      tk=twopi1/tu
      xk=twopi1/perl

      k=int(tdev/(dsi/v0))+1
      nthstep=k/(nlpoi-1)
      nx=nthstep*(nlpoi-1)
      nx=max(nx,nlpoi)

      if (mod(nx,nper).eq.0) then
        nlpoi=nlpoi-1
        nthstep=k/(nlpoi-1)
        nx=nthstep*(nlpoi-1)
        nx=max(nx,nlpoi)
      endif

      ds=xlen/dble(nx)
      dt=tdev/dble(nx)
      ds=dt*v0
      dx=xlen/dble(nx)

      yp=0.0d0

c      print*,nlpoi,nallo
      if (nlpoi.gt.nallo) then
        if (nallo.gt.0) then
          deallocate(tspl,ayrspl,ayispl,azrspl,azispl,w1,w2,w3,w4,w5,y2p,
     &      rarg,expo)
        endif
        allocate(amp(5,nlpoi),tspl(nlpoi),rarg(5,nlpoi),expo(2,nlpoi),
     &    ayrspl(nlpoi),ayispl(nlpoi),azrspl(nlpoi),azispl(nlpoi),
     &    w1(nlpoi),w2(nlpoi),w3(nlpoi),w4(nlpoi),w5(nlpoi),y2p(nlpoi))
        nallo=nlpoi
      endif

c      print*,'zmax:',za

      specnor=
     &  banwid
     &  /(4.0d0*pi1**2*clight1*hbarev1)
     &  /(4.0d0*pi1*eps01)
     &  *current/1.0d6 !per mm**2

      pownor=echarge1/16.0d0/pi1/pi1/eps01/clight1*current/1.0d6 !W/mm**2

      rspn=sqrt(specnor)

      if (idebug.ne.0) then
c        print*,"Calling urad_track"
      endif

      if (iobsv.eq.1) then
c        if (userchar(1).ne.'Ana') then
          call urad_track(
     &      icharge,
     &      gammai,dgamtot,
     &      xelec,yelec,zelec,vxelec,vyelec,vzelec,
     &      xf,yf,zf,efxn,efyn,efzn,
     &      xexit,yexit,zexit,vnxex,vnyex,vnzex,texit,ds,
     &      nstep,nlpoi,nthstep,traxyz,
     &      ieneloss,ivelofield,
     &      ifail,ith
     &      )
c        else
          dt=tdev/(nlpoi-1)
          dx=xlen/dble(nlpoi-1)
          betz0=park/gamma
          charge=dble(icharge)
c        endif
        if (nphener.gt.1) then
          dom=(phener(2)-phener(1))/hbarev1
        else
          dom=0.0d0
        endif
      endif !iobsv.eq.1

      ifr=0
      do ix=1,nlpoi

        if (userchar(1).ne.'Ana') then

          x=traxyz(1,ix)
          y=traxyz(2,ix)
          z=traxyz(3,ix)

          t=traxyz(4,ix)

          betx=traxyz(5,ix)*beta
          bety=traxyz(6,ix)*beta
          betz=traxyz(7,ix)*beta

        else

          dt=tdev/(nlpoi-1)
          t=dble(ix-1)*dt
c          t=traxyz(4,ix)
          st=sin(tk*t)
          ct=cos(tk*t)

          bety=0.0d0
          betz0=park/gamma*charge
          betz=betz0*st
          betx=beta*(1.0d0 - (betz/beta)**2/2.0d0 * (1.0d0+(betz/beta)**2/4.0d0))
          betx=beta*(1.0d0 - (betz/beta)**2/2.0d0)

          x=xelec+clight1*(
     &      2.0d0*ct*st*betz0**2*pi1*tu +
     &      4.0d0*beta**2*t-betz0**2*t)/4.0d0/beta
c          x=xelec+clight1*(
c     &      (4.0d0*ct*st**3*betz0**4*pi1*tu +
c     &      32.0d0*ct*st*beta**2*betz0**2*pi1*tu +
c     &      6.0d0*ct*st*betz0**4*pi1*tu +
c     &      64.0d0*beta**4*t-16.0d0*beta**2*betz0**2*t-3.0d0*betz0**4*t)/64.0d0/beta**3)
c          x=x*clight1+xelec
c          x=(x+xelec+t/tu*perl)/2.0d0
          y=0.0d0
          zamp=-betz0/tk*clight1
          z=zelec-zamp*ct
c          write(66,*)iobsv,ix,t,x,z,betx,betz,traxyz(1:4,ix),traxyz(5:7,ix)*beta

c          x=traxyz(1,ix)
c          y=traxyz(2,ix)
c          z=traxyz(3,ix)

c          t=traxyz(4,ix)

c          betx=traxyz(5,ix)*beta
c          bety=traxyz(6,ix)*beta
c          betz=traxyz(7,ix)*beta
        endif

        rnx=xobsv-x
        rny=yobsv-y
        rnz=zobsv-z

        h2=(rny*rny+rnz*rnz)/(rnx*rnx)

        if(h2.lt.0.01) then
          r=rnx*(1.0d0+(((((-0.0205078125D0*h2+0.02734375D0)*h2
     &      -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2)
        else
          r=sqrt(rnx**2+rny**2+rnz**2)
        endif

        rnx=rnx/r
        rny=rny/r
        rnz=rnz/r

        rarg(1,ix)=(rny*(rnx*bety-rny*betx)-rnz*(rnz*betx-rnx*betz))/r
        rarg(2,ix)=(rnz*(rny*betz-rnz*bety)-rnx*(rnx*bety-rny*betx))/r
        rarg(3,ix)=(rnx*(rnz*betx-rnx*betz)-rny*(rny*betz-rnz*bety))/r
        !phase0 noch unklar
        rarg(4,ix)=phase0+t-(rnx*x+rny*y+rnz*z)/clight1
        rarg(5,ix)=t

        tspl(ix)=t

      enddo !ix

      om=phener(1)/hbarev1-dom
      do ifr=1,nphener

        om=om+dom
        !om=phener(ifr)/hbarev1

        do ix=1,nlpoi

          rnnb(2:3)=rarg(2:3,ix)

          if (ifr.eq.1) then
            expo(1,ix)=cdexp(zi*om*rarg(4,ix))
            expo(2,ix)=cdexp(zi*dom*rarg(4,ix))
          else
            expo(1,ix)=expo(1,ix)*expo(2,ix)
          endif

          expom=expo(1,ix)
c          expom=cdexp(zi*om*rarg(4,ix))

          ayr=rnnb(2)*dreal(expom)
          ayi=rnnb(2)*dimag(expom)
          azr=rnnb(3)*dreal(expom)
          azi=rnnb(3)*dimag(expom)

          ayrspl(ix)=ayr
          ayispl(ix)=ayi
          azrspl(ix)=azr
          azispl(ix)=azi

          if (idebug.gt.1.and.yobsv.eq.0.0d0.and.zobsv.eq.0.0d0) then
            write(luno,*)ix,ifr,rarg(1:4,ix),om,tspl(ix),x,z,azr,azi
          endif

        enddo !nx

        eyr=0.0D0
        eyi=0.0D0
        ezr=0.0D0
        ezi=0.0D0

        if (user(7).eq.1) then

          do i=1,nlpoi-1
            dt2=(tspl(i+1)-tspl(i))/2.0d0
            eyr=eyr+dt2*(ayrspl(i)+ayrspl(i+1))
            eyi=eyi+dt2*(ayispl(i)+ayispl(i+1))
            ezr=ezr+dt2*(azrspl(i)+azrspl(i+1))
            ezi=ezi+dt2*(azispl(i)+azispl(i+1))
          enddo

        else if (user(7).eq.2) then

          yp1=9999.0d0
          ypn=9999.0d0

c          call util_spline_coef(tspl,ayrspl,nlpoi,yp1,ypn,y2p,w1,w2,w3,w4)
          do i=1,nlpoi-1
            dt=(tspl(i+1)-tspl(i))
            eyr=eyr
     &        +dt*0.5d0
     &        *(ayrspl(I)+ayrspl(I+1))
c     &        -dt**3/24.0d0
c     &        *(y2p(I)+y2p(I+1))
          enddo
c          call util_spline_coef(tspl,ayispl,nlpoi,yp1,ypn,y2p,w1,w2,w3,w4)
          do i=1,nlpoi-1
            dt=(tspl(i+1)-tspl(i))
            eyi=eyi
     &        +dt*0.5d0
     &        *(ayispl(I)+ayispl(I+1))
c     &        -dt**3/24.0d0
c     &        *(y2p(I)+y2p(I+1))
          enddo
c          call util_spline_coef(tspl,azrspl,nlpoi,yp1,ypn,y2p,w1,w2,w3,w4)
          do i=1,nlpoi-1
            dt=(tspl(i+1)-tspl(i))
            ezr=ezr
     &        +dt*0.5d0
     &        *(azrspl(I)+azrspl(I+1))
c     &        -dt**3/24.0d0
c     &        *(y2p(I)+y2p(I+1))
          enddo
c          call util_spline_coef(tspl,azispl,nlpoi,yp1,ypn,y2p,w1,w2,w3,w4)
          do i=1,nlpoi-1
            dt=(tspl(i+1)-tspl(i))
            ezi=ezi
     &        +dt*0.5d0
     &        *(azispl(I)+azispl(I+1))
c     &        -dt**3/24.0d0
c     &        *(y2p(I)+y2p(I+1))
          enddo

        else

          call util_spline_integral(tspl,ayrspl,nlpoi,eyr,y2p,w1,w2,w3,w4)
          call util_spline_integral(tspl,ayispl,nlpoi,eyi,y2p,w1,w2,w3,w4)
          call util_spline_integral(tspl,azrspl,nlpoi,ezr,y2p,w1,w2,w3,w4)
          call util_spline_integral(tspl,azispl,nlpoi,ezi,y2p,w1,w2,w3,w4)

        endif

        if (idebug.gt.1.and.yobsv.eq.0.0d0.and.zobsv.eq.0.0d0) then
          write(luna,*)ifr,phener(ifr),eyr,eyi,ezr,ezi,
     &      eyr**2+eyi**2+ezr**2+ezi**2
        endif

        aradey(ifr)=dcmplx((eyr),(eyi))*om
        aradez(ifr)=dcmplx((ezr),(ezi))*om

        enddo !nphener

      if (idebug.gt.1.and.yobsv.eq.0.0d0.and.zobsv.eq.0.0d0) then
        flush(luno)
        close(luno)
        flush(luna)
        close(luna)
      endif

      do kfreq=1,nphener

        aradex(kfreq)=aradex(kfreq)*rspn
        aradey(kfreq)=aradey(kfreq)*rspn
        aradez(kfreq)=aradez(kfreq)*rspn

        apolh=
     &    aradex(kfreq)*conjg(vstokes(1,1))
     &    +aradey(kfreq)*conjg(vstokes(1,2))
     &    +aradez(kfreq)*conjg(vstokes(1,3))

        apolr=
     &    aradex(kfreq)*conjg(vstokes(2,1))
     &    +aradey(kfreq)*conjg(vstokes(2,2))
     &    +aradez(kfreq)*conjg(vstokes(2,3))

        apoll=
     &    aradex(kfreq)*conjg(vstokes(3,1))
     &    +aradey(kfreq)*conjg(vstokes(3,2))
     &    +aradez(kfreq)*conjg(vstokes(3,3))

        apol45=
     &    aradex(kfreq)*conjg(vstokes(4,1))
     &    +aradey(kfreq)*conjg(vstokes(4,2))
     &    +aradez(kfreq)*conjg(vstokes(4,3))

        stok1=real(
     &    apolr*conjg(apolr)+
     &    apoll*conjg(apoll))

        stok2=-stok1+
     &    real(2.*apolh*conjg(apolh))

        stok3=
     &    real(2.*apol45*conjg(apol45))-
     &    stok1

        stok4=real(
     &    apolr*conjg(apolr)-
     &    apoll*conjg(apoll))

        stokes(1,kfreq)=stok1
        stokes(2,kfreq)=stok2
        stokes(3,kfreq)=stok3
        stokes(4,kfreq)=stok4

      enddo !nphener

      powden=powden*pownor

      if (idebug.ne.0) then
c        print*,"Leaving urad_spline_obspoi"
      endif
c+self.
      return
      end
