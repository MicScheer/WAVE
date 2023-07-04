*CMZ :  4.00/11 18/05/2021  19.13.31  by  Michael Scheer
*-- Author : Michael Scheer
C-----------------------------------------------------------
      subroutine bfourtabzy(xin,yin,zin,bxout,byout,bzout,axout,ayout,azout)

      use fbtabzymod

      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEND.

      integer :: nxbyfbtold,nxbzfbtold,nfourzyold,
     &  nxy,npow2,i,ifail,klinold,lunfft

      double precision
     &  xin,yin,zin,bxout,byout,bzout,axout,ayout,azout,
     &  bxf,byf,bzf,axf,ayf,azf,xiny,xinz

      complex, dimension(:), allocatable :: cffty,cfftz

      double precision, dimension(:), allocatable :: bypp,bzpp,
     &  ws1,ws2,ws3,ws4,bfft

      real, dimension(:), allocatable :: rbfft

      double precision xby1old,xbynold,by1old,bynold,
     &  xbz1old,xbznold,bz1old,bznold,dx,b,x,xi,xe,xc

      save

      if (inifbt.ne.0.or.nfourzyold.ne.nfourzy.or.
     &    nxbyfbt.ne.nxbyfbtold.or.xbyfbt(1).ne.xby1old.or.xbyfbt(nxbyfbt).ne.xbynold.or.
     &    nxbzfbt.ne.nxbzfbtold.or.xbzfbt(1).ne.xbz1old.or.xbzfbt(nxbzfbt).ne.xbznold.or.
     &    klinearfbt.ne.klinold
     &    ) then

        xi=xstart
        xc=xinter
        xe=xstop

        call fbtab_ini(xin,irbtab,irbtabzy,irbtabxyz)

        klinold=klinearfbt
        npow2=2**(nint(alog(float(nfourzy)/alog(2.)))+1)

        nfourzyold=nfourzy

        nxbyfbtold=nxbyfbt
        xby1old=xbyfbt(1)
        xbynold=xbyfbt(nxbyfbt)
        by1old=byfbt(1)
        bynold=byfbt(nxbyfbt)

        nxbzfbtold=nxbzfbt
        xbz1old=xbzfbt(1)
        xbznold=xbzfbt(nxbzfbt)
        bz1old=bzfbt(1)
        bznold=bzfbt(nxbzfbt)

        if (klinearfbt.eq.0) then

          if (iallofbt.ne.0) then
            deallocate(bypp,bzpp,ws1,ws2,ws3,ws4,cffty,cfftz,bfft,rbfft)
          endif

          nxy=max(nxbyfbt,nxbzfbt)
          allocate(bypp(nxbyfbt),bzpp(nxbzfbt),ws1(nxy),ws2(nxy),ws3(nxy),ws4(nxy))
          allocate(cffty(0:npow2/2-1),cfftz(0:npow2/2-1),
     &      bfft(npow2),rbfft(npow2))

          call util_spline_coef(xbyfbt,byfbt,nxbyfbt,0.0d0,0.0d0,bypp,
     &      ws1,ws2,ws3,ws4)

          dx=(xbyfbt(nxbyfbt)-xbyfbt(1))/dble(npow2)
          x=xbyfbt(1)
          call util_spline_inter(xbyfbt,byfbt,bypp,nxbyfbt,xbyfbt(1),bfft(1),-1)

          rbfft(1)=sngl(bfft(1))
          do i=2,npow2
            x=x+dx
            call util_spline_inter(xbyfbt,byfbt,bypp,nxbyfbt,x,bfft(i),0)
            rbfft(i)=sngl(bfft(i))
          enddo

          call util_rfft(npow2,rbfft,cffty,ifail)

          if (ifail.ne.0) then
            print*,"*** Error in bfourtabzy: Util_rfft failed ***"
            stop "*** Program WAVE aborted ***"
          endif

          open(newunit=lunfft,file=trim(filetb)//".fft")
          do i=0,npow2/2-1
            write(lunfft,*)i,real(cffty(i)),imag(cffty(i))
          enddo
          flush(lunfft)
          close(lunfft)

          call util_spline_coef(xbzfbt,bzfbt,nxbzfbt,0.0d0,0.0d0,bzpp,
     &      ws1,ws2,ws3,ws4)

          dx=(xbzfbt(nxbzfbt)-xbzfbt(1))/dble(npow2)
          x=xbzfbt(1)
          call util_spline_inter(xbzfbt,bzfbt,bzpp,nxbzfbt,xbzfbt(1),bfft(1),-1)

          rbfft(1)=sngl(bfft(1))
          do i=2,npow2
            x=x+dx
            call util_spline_inter(xbzfbt,bzfbt,bzpp,nxbzfbt,x,bfft(i),0)
            rbfft(i)=sngl(bfft(i))
          enddo

          call util_rfft(npow2,rbfft,cfftz,ifail)

          if (ifail.ne.0) then
            print*,"*** Error in bfourtabzy: Util_rfft failed ***"
            stop "*** Program WAVE aborted ***"
          endif

          open(newunit=lunfft,file=trim(filetbz)//".fft")
          do i=0,npow2/2-1
            write(lunfft,*)i,real(cfftz(i)),imag(cfftz(i))
          enddo
          flush(lunfft)
          close(lunfft)

          iallofbt=1

        else !klinearfbt

          if (iallofbt.ne.0) then
            deallocate(bypp,bzpp,ws1,ws2,ws3,ws4,cffty,cfftz,bfft,rbfft)
          endif

          nxy=max(nxbyfbt,nxbzfbt)
          allocate(bypp(nxbyfbt),bzpp(nxbzfbt),ws1(nxy),ws2(nxy),ws3(nxy),ws4(nxy))
          allocate(cffty(0:npow2/2-1),cfftz(0:npow2/2-1),
     &      bfft(npow2),rbfft(npow2))

          dx=(xbyfbt(nxbyfbt)-xbyfbt(1))/dble(npow2)
          x=xbyfbt(1)-dx
          do i=1,npow2
            x=x+dx
            call util_interpol_linear(nxbyfbt,xbyfbt,byfbt,x,bfft(i),ifail)
            if (ifail.ne.0) then
              print*,"*** Error in bfourtabzy: Util_interpol_linear failed ***"
              stop "*** Program WAVE aborted ***"
            endif
            rbfft(i)=sngl(bfft(i))
          enddo

          call util_rfft(npow2,rbfft,cffty,ifail)

          dx=(xbyfbt(nxbzfbt)-xbzfbt(1))/dble(npow2)
          x=xbzfbt(1)-dx
          do i=1,npow2
            x=x+dx
            call util_interpol_linear(nxbzfbt,xbzfbt,bzfbt,x,bfft(i),ifail)
            if (ifail.ne.0) then
              print*,"*** Error in bfourtabzy: Util_interpol_linear failed ***"
              stop "*** Program WAVE aborted ***"
          endif
            rbfft(i)=sngl(bfft(i))
          enddo

          call util_rfft(npow2,rbfft,cfftz,ifail)

          if (ifail.ne.0) then
            print*,"*** Error in bfourtabzy: Util_rfft failed ***"
            stop "*** Program WAVE aborted ***"
          endif

          if (ifourzy0.eq.11) then
            cffty(0)=(0.0,0.0)
            cfftz(0)=(0.0,0.0)
          else if (ifourzy0.eq.10) then
            cfftz(0)=(0.0,0.0)
          else if (ifourzy0.eq.1) then
            cffty(0)=(0.0,0.0)
          endif

          iallofbt=1

        endif !klinearfbt

        if (xi.eq.9999.0d0) xstart=dmin1(xbyfbt(1),xbzfbt(1))
        if (xe.eq.9999.0d0) xstop=dmax1(xbyfbt(nxbyfbt),xbzfbt(nxbzfbt))

      endif !inifbt

      bxout=0.0d0
      byout=0.0d0
      bzout=0.0d0

      axout=0.0d0
      ayout=0.0d0
      azout=0.0d0

      xiny=xin+fourysh

      if (xiny.ge.xbyfbt(1).and.xiny.le.xbyfbt(nxbyfbt)) then
        call bfoursincos(xiny-xbyfbt(1),yin,zin,
     &    bxout,byout,bzout,axout,ayout,azout,
     &    npow2/2,xbyfbt(nxbyfbt)-xbyfbt(1),0.0d0,cffty,inifbt)
      endif

      xinz=xin+fourzsh

      if (xinz.ge.xbzfbt(1).and.xinz.le.xbzfbt(nxbzfbt)) then
        call bfoursincos(xinz-xbzfbt(1),yin,zin,
     &    bxf,byf,bzf,axf,ayf,azf,
     &    npow2/2,xbzfbt(nxbzfbt)-xbzfbt(1),0.0d0,cfftz,0)
        bxout=bxout+bxf
        byout=byout+bzf
        bzout=bzout+byf
        axout=bxout+axf
        ayout=byout+azf
        azout=bzout+ayf
      endif

      return
      end
