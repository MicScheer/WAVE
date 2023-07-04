*CMZ : 00.00/15 24/10/2012  14.36.19  by  Michael Scheer
*CMZ : 00.00/11 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/07 02/05/2008  13.10.35  by  Michael Scheer
*CMZ : 00.00/06 01/11/2007  18.40.58  by  Michael Scheer
*-- Author :    Michael Scheer   31/10/2007
      subroutine util_rfft_spline(npoi,nfour,x,y,a0,ac,as,mode,ifail)

c +PATCH,//UTIL/FOR
c +DECK,util_rfft_spline.

c Input and Output in double precision, but fft is only in single precision!

      implicit none

      integer ndimp
      parameter (ndimp=2**15)

      double precision x(npoi),y(npoi),a0,ac(npoi),as(npoi),dx
      integer npoi,ifail,npow,nfour,i,nfouro,mode

      double precision
     &  ws1(ndimp),ws2(ndimp),ws3(ndimp),ws4(ndimp),ws5(ndimp),
     &  xs(ndimp),ys(ndimp),y2p(ndimp)

c Neu, weil ifort auf dinux die Uebergabe von rfft nicht akzeptiert.
c WENN DAS MAL GUTGEHT...
      complex cfft(0:ndimp/2-1)
      real yr(ndimp),rfft(ndimp)
      equivalence(cfft,rfft)

      if (npoi.lt.2.or.npoi.gt.ndimp) then
        ifail=-1
        return
      endif

      if (nfour.lt.2) then
        ifail=1
        nfour=npoi
      endif

      ac=0.0d0
      as=0.0d0

      npow=int(alog(float(npoi))/alog(2.0))
      if (2**(npow+1).le.npoi) npow=npow+1

      nfouro=nfour
      nfour=2**nint(alog(float(nfouro))/alog(2.0))
      nfour=min(nfour,2**npow)

      if (nfour.gt.nfouro) then
        nfour=nfour/2
      endif

      call util_spline_coef_periodic(x,y,npoi,y2p,ws1,ws2,ws3,ws4,ws5,ifail)
      if (ifail.ne.0) return

      dx=(x(npoi)-x(1))/(nfour-1)
      xs(1)=x(1)
      call util_spline_inter(x,y,y2p,npoi,xs(1),ys(1),-1)

      do i=2,nfour-1
        xs(i)=xs(i-1)+dx
        call util_spline_inter(x,y,y2p,npoi,xs(i),ys(i),0)
      enddo !nfour

      xs(nfour)=x(npoi)
      call util_spline_inter(x,y,y2p,npoi,xs(nfour),ys(nfour),0)

      if (mode.ne.0) then
c shuffle data such that x is centered, i.e. phase-shift by pi
        yr(1:nfour/2)=ys(nfour/2+1:nfour)
        yr(nfour/2+1:nfour)=ys(1:nfour/2)
      else !mode
        yr=ys
      endif !mode

      call util_rfft(nfour,yr,cfft,ifail)
      if (ifail.ne.0) return

      a0=rfft(1)

      do i=1,nfour/2-1
        ac(i)=rfft(2*i+1)
        as(i)=rfft(2*i+2)
      enddo

      return
      end
