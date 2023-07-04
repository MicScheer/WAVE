*CMZ : 00.00/15 26/06/2012  15.42.50  by  Michael Scheer
*CMZ : 00.00/14 27/09/2011  12.35.34  by  Michael Scheer
*CMZ : 00.00/13 09/09/2011  09.25.17  by  Michael Scheer
*CMZ : 00.00/12 01/09/2011  12.02.02  by  Michael Scheer
*CMZ : 00.00/11 02/08/2011  08.49.19  by  Michael Scheer
*CMZ : 00.00/10 29/12/2010  20.54.24  by  Michael Scheer
*CMZ : 00.00/07 19/08/2010  12.20.18  by  Michael Scheer
*CMZ : 00.00/06 17/10/2007  16.49.18  by  Michael Scheer
*CMZ : 00.00/05 25/12/2006  14.48.59  by  Michael Scheer
*CMZ : 00.00/04 07/12/2006  16.26.31  by  Michael Scheer
*CMZ : 00.00/03 24/11/2006  11.13.44  by  Michael Scheer
*CMZ : 00.00/02 21/11/2006  15.26.03  by  Michael Scheer
      subroutine util_smooth(mode,mord,nival,kval,lunxfit,ifail,cerr)

c *** Die Wichtung ist fraglich, da zumindest im Clustermode die Gewichte nur
c innerhalb des Intervalles wirken!?


c Input arguments (see also util_smooth.nam for util_smooth_main.f):
c
c      mode:   mode=1 !1: fit polynom
c            else: cluster mode,
c            i.e. xcen(ival),ycen(ival) are center of gravity (mord=1)
c
c      kval:   1: get xfit-values from filexfit (common-block)
c             -1: under control of MINUIT
c            2: some number of points per interval (within +/- nival-1)
c               number of interval is given is nival
c            3: equal intervals
c
c      nival:  number of intervals for kval.ne.1
c
c      if ifixpar.ne.0, fixed parameters are read from util_smooth.fix
c
c      lunxfit: LUN for file filexfit
c
c      ifail:  Status flag
c
c      cerr (character*128): Error code for ifail.ne.0
c
c

c    Input in common/usmooth/
c    ndata,xdata,ydata,edata,filexfit
c

c Output:
c
c      xfit, yfit, vfit. yfit are the fitted y-values yfit(xfit),
c      while vfit contains all fitted parameters
c
c      ysmooth(xdata),xcen,ycen,npar,nival
c      ysmooth is only calculated within region of fit
c
c      chi2 or chi2/ndf
c

c Comment:
c
c     edata are the errors (1 sigma) of ydata, i.e. ydata = ydata +/- edata
c     chi2 = sum( (f(xdata(i))-ydata(i)/edata(i))**2 ,
c     i = 1, all data within [xfit(1),xfit(nival)]),
c     i.e. all contributions to the equation-system are weighted
c     by 1/edata(i)**2.
c
c     In the case of npar=ndata, the cubic-spline-fit does not agree exactliy
c     to the natural spline, since also the second derivatives are fitted for
c     the first and last data point.
c
c
c     Interval bounderies should coincide with xdata, otherwise strong
c     oscillations may occure for interpolated data. kval=2 is recommended.
c     On the other hand, if boundaries and xdata coninced, the determinate of
c     the equation system may become zero and has only a solution due to
c      numerical errors. So, it could be a good idea, to
c     test values with an epsilon near xdata.
c
c
      implicit none

*KEEP,UTIL_SMOOTH.
      integer ndatap,nfitp,nivalp

      parameter(ndatap=100000,nfitp=100,nivalp=100)

      integer nvarp,nvar
      parameter (nvarp=100)

      double precision
     &  fchi2,fchi2min,
     &  fvarfin(nvarp),fvar(nvarp),fparopt(nvarp),fweight(nvarp)

      double precision
     &  chi2,chi2ndf,
     &  xfit(nivalp+1),yfit(nivalp+1),
     &  vfit(nivalp+1),
     &  xcen(nivalp),ycen(nivalp),
     &  xdata(ndatap),ydata(ndatap),edata(ndatap),
     &  ysmooth(ndatap),xexe,xexa,xlow,xhig,ylow,yhig,
     &  coef(nivalp+1),dspps1,xknots(nivalp+1),
     &  ws(ndatap*(nivalp+5)+nivalp*(nivalp+1))

      integer ndata,nfit,mdatival,nexpand,
     &  interval(nivalp+1),ninterval(nivalp+1),ifixpar
     &  ,modef,kvalf,mordf,nivalf,ifailf,lunxfitf,
     &  nknots,nw,knotmode

      character(256) filexfit
      character(128) cerrf

      common/usmooth/
     &  chi2,chi2ndf,
     &  xfit,yfit,vfit,
     &  xcen,ycen,
     &  xdata,ydata,edata,
     &  ysmooth,xexa,xexe,xlow,xhig,ylow,yhig,
     &  ndata,nfit,mdatival,interval,ninterval,nexpand,ifixpar,
     &  filexfit,
     &  fchi2,fchi2min,fvarfin,fvar,fparopt,fweight,nvar
     &  ,modef,kvalf,mordf,nivalf,ifailf,lunxfitf,cerrf,
     &  nknots,coef,xknots,nw,ws,knotmode
*KEND.

      double precision a,b,
     &  ab(nivalp+1),a2(nivalp+1),b2(nivalp+1),ay(nivalp),by(nivalp),
     &  fmat(9,nivalp+1),
     &  fmat6(6,nivalp+1),
     &  fmat7(7,nivalp+1),
     &  abi,b2i,abm,aym,byi,a2m,
     &  dxfit,dx,y,
     &  xfiti,x,dx2,
     &  xi2,xi1,
     &  cay(nivalp),caaay(nivalp),aaa,bbb,h6,
     &  cby(nivalp),cbbby(nivalp),
     &  ca2(nivalp),cb2(nivalp),cab(nivalp),
     &  caaaa(nivalp),cbbbb(nivalp),cabbb(nivalp),cbaaa(nivalp),
     &  caaabbb(nivalp),caaa2(nivalp),cbbb2(nivalp),errrms,
     &  weight,wsum,dxi,dxi2,dxo,dxo2,cy,cyp,cypp,
     &  cyy(nivalp),cyyp(nivalp),cyypp(nivalp),cypyp(nivalp),cypypp(nivalp),
     &  cyppypp(nivalp),cyd(nivalp),cypd(nivalp),cyppd(nivalp),
     &  fixpar(nivalp)

      integer mode,mord,kval,ifail,i,nival,ival,ival1,ivalm,ivalp,lunxfit,
     &  idata,ndata1,ndata2,inter1,inter2,npar,ipar,i1,i2,i3,i4,i5,i6,
     &  i7,i8,l1,l2,nord,mdata,kminuit

      character(128) cerr

      data kminuit/0/

      if (kval.lt.0) then
        kminuit=kval
        kval=1
      endif

      ifail=0
      cerr='--- util_smooth returned ok'

      if (mord.lt.1) then
        cerr='*** Error in util_smooth: Bad order'
        ifail=1
        return
      endif

      if (mode.eq.-4) then
        kminuit=0
        mode=4
      endif

      if (kminuit.eq.0) then

        if (mord.eq.1) then
          npar=nival+1
          nord=1
        else if (mord.eq.2) then
          npar=nival+2
          nord=2
        else if (mord.eq.3) then
          npar=2*nival+2
          nord=2
        else if (mord.eq.4) then

          if (kval.eq.1) then

            nival=0
            open(unit=lunxfit,file=filexfit)
111         nival=nival+1
            read(lunxfit,*,end=999)xfit(nival)
            goto 111
999         nival=nival-2
            close(lunxfit)
            call util_sort(nival+1,xfit)

            if (mord.eq.4) then
              knotmode=3
              nknots=nival+1
              xknots(1:nknots)=xfit(1:nknots)
            endif

          endif !kval=1

          nord=1

          goto 123

        else !mord
          cerr='*** Error in util_smooth: Bad order'
          ifail=1
          return
        endif

        fixpar=-9999.0d0
        if (ifixpar.ne.0) then
          open(unit=99,file='util_smooth.fix')
7         continue
          call util_skip_comment_end(99,ifixpar)
          if (ifixpar.eq.1) goto 77
          read(99,*,end=77)ifixpar,fixpar(ifixpar)
          goto 7
77        close(99)
        endif

      else !kminuit

        if (mord.eq.4) then
          xknots(1:nknots)=xfit(1:nknots)
          goto 123
        endif

      endif !kminuit

c{ define smoothing intervals

      if (kval.eq.1) then

        if (kminuit.eq.0) then

          nival=0

          open(unit=lunxfit,file=filexfit)

11        nival=nival+1
          read(lunxfit,*,end=99)xfit(nival)
          goto 11
99        nival=nival-2
          close(lunxfit)

          call util_sort(nival+1,xfit)

          if (mord.eq.4) then
            knotmode=3
            nknots=nival+1
            xknots(1:nknots)=xfit(1:nknots)
            kval=1
          endif

        endif !kminuit

        if (nival.lt.1) then
          cerr='*** Error in util_smooth: Bad number of intervals'
          ifail=1
          return
        endif

        if (nival.gt.nivalp) then
          cerr='*** Error in util_smooth: Dimension of intervals exceeded'
          ifail=1
          return
        endif

        call util_sort(nival+1,xfit)

        if (mord.eq.1) then
          npar=nival+1
          nord=1
        else if (mord.eq.2) then
          npar=nival+2
          nord=2
        else if (mord.eq.3) then
          npar=2*nival+2
          nord=2
        else if (mord.eq.4) then
          nknots=nival+1
          if (kval.eq.3) nknots=max(6,nknots)
          npar=nknots
          nord=1
        else !mord
          cerr='*** Error in util_smooth: Bad order'
          ifail=1
          return
        endif

        ninterval(nival+1)=0

        ndata1=0
        do idata=1,ndata-1
          if (xdata(idata).ge.xfit(1).and.ndata1.eq.0) then
            ndata1=idata
          endif
          if (xdata(idata).le.xfit(nival+1)) then
            ndata2=idata
          endif
          if (xdata(idata).ge.xdata(idata+1)) then
            cerr=
     & '*** Error in util_smooth: Bad x-data, must be strictly raising'
            ifail=1
            return
          endif
        enddo

        if (xdata(ndata).le.xfit(nival+1)) then
          ndata2=ndata
        endif

        idata=ndata1
        do ival=1,nival-1
12        if (xdata(idata).lt.xfit(ival+1)) then
            ninterval(ival)=ninterval(ival)+1
            if (interval(ival).eq.0) interval(ival)=idata
            idata=idata+1
            if (idata.gt.ndata2) then
              goto 92
            else
              goto 12
            endif
          endif
        enddo !ival

        ival=nival
        if (idata.eq.0) then
          cerr='*** Error in util_smooth: Too few data in interval:'
          ifail=ival
          return
        endif

121    if (xdata(idata).le.xfit(ival+1)) then
          ninterval(ival)=ninterval(ival)+1
          if (interval(ival).eq.0) interval(ival)=idata
          idata=idata+1
          if (idata.le.ndata2) goto 121
        endif

92      interval(nival+1)=interval(nival)+ninterval(nival)

c data per interval

        if (mord.eq.1) then
          npar=nival+1
          nord=1
          do ival=1,nival
            if (ninterval(ival).lt.nord) then
              cerr='*** Error in util_smooth: Too few data in interval:'
              ifail=ival
              return
            endif
          enddo
        else if (mord.eq.2) then
          npar=nival+2
          nord=2
          do ival=1,nival
            if (ninterval(ival).lt.1) then
              cerr='*** Error in util_smooth: Too few data in interval:'
              ifail=ival
              return
            endif
          enddo
        else if (mord.eq.3) then
          npar=2*nival+2
          nord=2
          do ival=1,nival
            if (ninterval(ival).lt.nord) then
              cerr='*** Error in util_smooth: Too few data in interval:'
              ifail=ival
              return
            endif
          enddo
        else if (mord.eq.4) then
          nknots=nival+1
          npar=nknots
          nord=1

        else !mord
          cerr='*** Error in util_smooth: Bad order'
          ifail=1
          return
        endif

      else if (kval.eq.2) then

        if (nival.lt.1) then
          cerr='*** Error in util_smooth: Bad number of intervals'
          ifail=1
          return
        endif

c data per interval

        mdatival=ndata/nival

        if (mdatival.lt.mord) then
          cerr='*** Error in util_smooth: Too few data per interval'
          ifail=1
          return
        endif

c simply same number of data per interval

        xfit(1)=xdata(1)
        xfit(nival+1)=xdata(ndata)
        ninterval(nival)=ndata-mdatival*(nival-1)
        interval(1)=1
        interval(nival+1)=ndata+1

        do ival=1,nival-1
          i=ival*mdatival
          ninterval(ival)=mdatival
          interval(ival+1)=interval(ival)+mdatival
          xfit(ival+1)=(xdata(i)+xdata(i+1))/2.0d0
        enddo

        ndata1=1
        ndata2=ndata

      else if (kval.eq.3) then

        if (nival.lt.1) then
          cerr='*** Error in util_smooth: Bad number of intervals'
          ifail=1
          return
        endif

        if (nival.gt.nivalp) then
          cerr='*** Error in util_smooth: Dimension of intervals exceeded'
          ifail=1
          return
        endif

        dxfit=(xdata(ndata)-xdata(1))/nival
        xfit(1)=xdata(1)
        xfit(nival+1)=xdata(ndata)

        do ival=1,nival
          interval(ival)=0
          ninterval(ival)=0
          xfit(ival)=xfit(1)+(ival-1)*dxfit
        enddo

        ninterval(nival+1)=0

        ndata1=0
        do idata=1,ndata-1
          if (xdata(idata).ge.xfit(1).and.ndata1.eq.0) then
            ndata1=idata
          endif
          if (xdata(idata).le.xfit(nival+1)) then
            ndata2=idata
          endif
          if (xdata(idata).ge.xdata(idata+1)) then
            cerr=
     & '*** Error in util_smooth: Bad x-data, must be strictly raising'
            ifail=1
            return
          endif
        enddo

        if (xdata(ndata).le.xfit(nival+1)) then
          ndata2=ndata
        endif

        idata=ndata1
        do ival=1,nival-1
122       if (xdata(idata).lt.xfit(ival+1)) then
            ninterval(ival)=ninterval(ival)+1
            if (interval(ival).eq.0) interval(ival)=idata
            idata=idata+1
            if (idata.gt.ndata2) then
              goto 922
            else
              goto 122
            endif
          endif
        enddo !ival

        ival=nival
1212    if (xdata(idata).le.xfit(ival+1)) then
          ninterval(ival)=ninterval(ival)+1
          if (interval(ival).eq.0) interval(ival)=idata
          idata=idata+1
          if (idata.le.ndata2) goto 1212
        endif

922     interval(nival+1)=interval(nival)+ninterval(nival)

c data per interval

        do ival=1,nival
          if (ninterval(ival).lt.nord) then
            cerr='*** Error in util_smooth: Too few data in interval:'
            ifail=ival
            return
          endif
        enddo

      else !kval

        if (kminuit.eq.0) then
          cerr='*** Error in util_smooth: Bad kval'
          ifail=1
          return
        endif

      endif !kval

c} define smoothing intervals

      if (kminuit.ne.0) then
        goto 123
      endif

      mdata=ndata2-ndata1+1

      errrms=0.0d0
      do idata=ndata1,ndata2
        if (ydata(idata).ne.0.0d0) then
          errrms=errrms+(edata(idata)/ydata(idata))**2
        else
          errrms=errrms+(edata(idata))**2
        endif
      enddo

      errrms=sqrt(errrms/mdata)

      if (errrms.ne.0.0d0) then
        do idata=ndata1,ndata2
          if (edata(idata).eq.0.0d0) then
            edata(idata)=ydata(idata)*errrms
          endif
        enddo
        errrms=0.0d0
        do idata=1,ndata
          errrms=errrms+(edata(idata))**2
      enddo
      errrms=sqrt(errrms/mdata)
      else
        do idata=1,ndata
          if (edata(idata).eq.0.0d0) then
            edata(idata)=1.0d0
          endif
        enddo
      endif

      if (nival*nord+1.gt.nfitp) then
          cerr='*** Error in util_smooth: Dimension of fit parameters exceeded'
          ifail=1
          return
        endif

      if (mdata.lt.npar) then
        cerr='*** Error in util_smooth: Less data than parameters to fit:'
        ifail=ival
        return
      endif

123   continue

      if (mord.eq.1) then

        if (mode.eq.1) then

c{ first order, i.e. polygon

          do ival=1,nival

            ival1=ival+1

            a2(ival)=0.0d0
            b2(ival)=0.0d0
            ab(ival)=0.0d0
            ay(ival)=0.0d0
            by(ival)=0.0d0

            dxfit=xfit(ival1)-xfit(ival)

            do i=interval(ival),interval(ival1)-1

              weight=1.0d0/edata(i)

              a=(xdata(i)-xfit(ival))/dxfit*weight
              b=weight-a
              y=ydata(i)*weight

              a2(ival)=a2(ival)+a*a
              b2(ival)=b2(ival)+b**2
              ab(ival)=ab(ival)+b*a
              ay(ival)=ay(ival)+a*y
              by(ival)=by(ival)+b*y

            enddo

          enddo !ival

          do ival=1,nival+1

            ivalm=ival-1
            ivalp=ival+1

            if (ival.ge.1.and.ival.le.nival) then
              b2i=b2(ival)
              abi=ab(ival)
              byi=by(ival)
            else
              b2i=0.0d0
              abi=0.0d0
              byi=0.0d0
            endif

            if (ivalm.ge.1.and.ivalm.le.nival) then
              abm=ab(ivalm)
              a2m=a2(ivalm)
              aym=ay(ivalm)
            else
              abm=0.0d0
              a2m=0.0d0
              aym=0.0d0
            endif

c non-vanishing matrix elements

            fmat(1,ival)=abm
            fmat(2,ival)=b2i+a2m
            fmat(3,ival)=abi

c inhomogenity

            fmat(4,ival)=byi+aym

          enddo !ival

c solution

          ival=1
          fmat(3,ival)=fmat(3,ival)/fmat(2,ival)
          fmat(4,ival)=fmat(4,ival)/fmat(2,ival)
          fmat(2,ival)=1.0d0

          do ival=2,nival+1

            ivalm=ival-1

            if (fmat(1,ival).ne.0.0d0) then

              fmat(2,ival)=fmat(2,ival)/fmat(1,ival)-fmat(3,ivalm)
              fmat(3,ival)=fmat(3,ival)/fmat(1,ival)
              fmat(4,ival)=fmat(4,ival)/fmat(1,ival)-fmat(4,ivalm)
              fmat(1,ival)=0.0d0

              fmat(3,ival)=fmat(3,ival)/fmat(2,ival)
              fmat(4,ival)=fmat(4,ival)/fmat(2,ival)
              fmat(2,ival)=1.0d0

            endif

          enddo !ival

          do ival=nival,1,-1

            ivalp=ival+1

            if (fmat(3,ival).ne.0.0d0) then

              fmat(2,ival)=fmat(2,ival)/fmat(3,ival)
              fmat(4,ival)=fmat(4,ival)/fmat(3,ival)-fmat(4,ivalp)
              fmat(3,ival)=0.0d0

              fmat(4,ival)=fmat(4,ival)/fmat(2,ival)
              fmat(2,ival)=1.0d0

            endif

            yfit(ival)=fmat(4,ival)

          enddo !ival

          yfit(nival+1)=fmat(4,nival+1)

          do ival=1,nival+1
            vfit(ival)=yfit(ival)
          enddo

c{ smoothed data

          do ival=1,nival

            ivalp=ival+1

            a=(yfit(ivalp)-yfit(ival))/(xfit(ivalp)-xfit(ival))

            do i=interval(ival),interval(ivalp)-1
              ysmooth(i)=yfit(ival)+a*(xdata(i)-xfit(ival))
            enddo

            xcen(ival)=(xfit(ivalp)+xfit(ival))/2.0d0
            ycen(ival)=yfit(ival)+a*(xcen(ival)-xfit(ival))

          enddo !ival

        else ! (mode.eq.1) then

c{ cluster mode

          if (nival.lt.2) then
            cerr=
     &        '*** Error in util_smooth: Cluster mode, but only less then 2 intervals'
            ifail=1
            return
          endif

          do ival=1,nival

            ival1=ival+1

            xcen(ival)=0.0d0
            ycen(ival)=0.0d0

            wsum=0.0d0

            do i=interval(ival),interval(ival1)-1
              weight=1.0d0/edata(i)
              xcen(ival)=xcen(ival)+xdata(i)*weight
              ycen(ival)=ycen(ival)+ydata(i)*weight
              wsum=wsum+weight
            enddo !idata

            xcen(ival)=xcen(ival)/wsum
            ycen(ival)=ycen(ival)/wsum

          enddo !ival

          yfit(1)=ycen(1)+
     &      (ycen(2)-ycen(1))/(xcen(2)-xcen(1))*
     &      (xfit(1)-xcen(1))
          yfit(nival+1)=ycen(nival)+
     &      (ycen(nival)-ycen(nival-1))/(xcen(nival)-xcen(nival-1))*
     &      (xfit(nival+1)-xcen(nival))

          do i=2,nival
              ivalm=i-1
              ivalp=i
              yfit(i)=ycen(ivalm)+
     &          (ycen(ivalp)-ycen(ivalm))/(xcen(ivalp)-xcen(ivalm))*
     &          (xfit(i)-xcen(ivalm))
          enddo !ival

          do ival=1,nival

            ivalp=ival+1

            a=(yfit(ivalp)-yfit(ival))/(xfit(ivalp)-xfit(ival))

            do i=interval(ival),interval(ivalp)-1
              ysmooth(i)=yfit(ival)+a*(xdata(i)-xfit(ival))
            enddo

          enddo !ival

c} cluster mode

        endif !(mode.eq.1) then

c} smoothed data

c} first order, i.e. polygon

      else if (mord.eq.2) then

c{ second order
c function f1=sum(i=1:ndata_1; a1i + b1 * (x(i)-xi1) + c1 (x(i)-xi1)**2)

        if (ifixpar.ne.0) then
          cerr='*** Warning: ifixpar not yet implemented for second order spline ***'
          ifail=1
          return
        endif

        if (nival*3.gt.nivalp) then
          cerr='*** Error in util_smooth: Dimension of intervals exceeded'
          ifail=1
          return
        endif

        cyy(1:nival)=0.0d0
        cyyp(1:nival)=0.0d0
        cyypp(1:nival)=0.0d0

        cypyp(1:nival)=0.0d0
        cypypp(1:nival)=0.0d0

        cyppypp(1:nival)=0.0d0

        cyd(1:nival)=0.0d0
        cypd(1:nival)=0.0d0
        cyppd(1:nival)=0.0d0

        ival=1

        i1=1+(ival-1)*3
        i2=i1+1

        inter1=interval(i1)
        inter2=interval(i2)-1

        xi1=xfit(i1)
        xi2=xfit(i2)

        dxi=xi2-xi1
        dxi2=dxi*dxi

        do i=inter1,inter2

          weight=1.0d0/edata(i)

          x=xdata(i)
          y=ydata(i)

          dx=x-xi1
          dx2=dx*dx

          cy=1.0d0
          cyp=dx
          cypp=dx2

          cyy(ival)   =cyy(ival) + cy*cy*weight
          cyyp(ival)  =cyyp(ival) + cy*cyp*weight
          cyypp(ival) =cyypp(ival)+ cy*cypp*weight

          cypyp(ival)  =cypyp(ival) + cyp*cyp*weight
          cypypp(ival) =cypypp(ival)+ cyp*cypp*weight

          cyppypp(ival) =cyppypp(ival)+ cypp*cypp*weight

          cyd(ival)=cyd(ival)+y*cy*weight
          cypd(ival)=cypd(ival)+y*cyp*weight
          cyppd(ival)=cyppd(ival)+y*cypp*weight

        enddo  !i (data)

        fmat6(1:6,1:nival*3)=0.0d0

        i1=ival
        i2=i1+1
        i3=i2+1

        fmat6(3,i1)=cyy(ival)
        fmat6(4,i1)=cyyp(ival)
        fmat6(5,i1)=cyypp(ival)

        fmat6(2,i2)=cyyp(ival)
        fmat6(3,i2)=cypyp(ival)
        fmat6(4,i2)=cypypp(ival)

        fmat6(1,i3)=cyypp(ival)
        fmat6(2,i3)=cypypp(ival)
        fmat6(3,i3)=cyppypp(ival)

        vfit(i1)=cyd(ival)
        vfit(i2)=cypd(ival)
        vfit(i3)=cyppd(ival)

        dxo=dxi
        dxo2=dxi2

        do ival=2,nival

          i3=(ival-1)*3
          i4=i3+1
          i5=i4+1
          i6=i5+1

          l1=ival
          l2=ival+1

          inter1=interval(l1)
          inter2=interval(l2)-1

          xi1=xfit(l1)
          xi2=xfit(l2)

          dxi=xi2-xi1
          dxi2=dxi*dxi

          do i=inter1,inter2

            weight=1.0d0/edata(i)

            x=xdata(i)
            y=ydata(i)

            dx=x-xi1
            dx2=dx*dx

            cy=1.0d0
            cyp=dx
            cypp=dx2

            cyy(ival)   =cyy(ival) + cy*cy*weight
            cyyp(ival)  =cyyp(ival) + cy*cyp*weight
            cyypp(ival) =cyypp(ival)+ cy*cypp*weight

            cypyp(ival)  =cypyp(ival) + cyp*cyp*weight
            cypypp(ival) =cypypp(ival)+ cyp*cypp*weight

            cyppypp(ival) =cyppypp(ival)+ cypp*cypp*weight

            cyd(ival)=cyd(ival)+y*cy*weight
            cypd(ival)=cypd(ival)+y*cyp*weight
            cyppd(ival)=cyppd(ival)+y*cypp*weight

          enddo  !i (data)

          fmat6(1,i3)=1.0d0
          fmat6(2,i3)=dxo
          fmat6(3,i3)=dxo2
          fmat6(4,i3)=-1.0d0
          vfit(i3)=0.0d0

          fmat6(1,i4)=1.0d0
          fmat6(2,i4)=2.0d0*dxo
          fmat6(4,i4)=-1.0d0
          vfit(i4)=0.0d0

          fmat6(2,i5)=cyyp(ival)
          fmat6(3,i5)=cypyp(ival)
          fmat6(4,i5)=cypypp(ival)
          vfit(i5)=cypd(ival)

          fmat6(1,i6)=cyypp(ival)
          fmat6(2,i6)=cypypp(ival)
          fmat6(3,i6)=cyppypp(ival)
          vfit(i6)=cyppd(ival)

        enddo !ival

        call util_smooth_gauss_6x6_s(nival*3,fmat6,vfit)

c} second order

c{ smoothed data

        do ival=1,nival

          xi1=xfit(ival)
          ipar=3*(ival-1)+1

          do i=interval(ival),interval(ival+1)-1

            x=xdata(i)

            dx=x-xi1
            dx2=dx**2

            ysmooth(i)=vfit(ipar)+vfit(ipar+1)*dx+vfit(ipar+2)*dx**2

          enddo

          xcen(ival)=(xfit(ival)+xfit(ival+1))/2.0d0
          x=xcen(ival)

          dx=x-xi1
          dx2=dx**2

          ycen(ival)=vfit(ipar)+vfit(ipar+1)*dx+vfit(ipar+2)*dx2
          yfit(ival)=vfit(ipar)

        enddo !ival

        x=xfit(nival+1)

        dx=x-xi1
        dx2=dx**2

        yfit(nival+1)=
     &    vfit(ipar)+vfit(ipar+1)*dx+vfit(ipar+2)*dx2

      else if (mord.eq.3) then

c{ third order, i.e. cubic-spline

c function values and second derivatives are fitted on intervall boundaries

        if (nival*4.gt.nivalp) then
          cerr='*** Error in util_smooth: Dimension of intervals exceeded'
          ifail=1
          return
        endif

        if (ifixpar.ne.0) then
          cerr='*** Warning: ifixpar not yet implemented for cubic spline ***'
          ifail=1
          return
        endif

        do ival=1,nival

          ival1=ival+1

          ca2(ival)=0.0d0
          cab(ival)=0.0d0
          cb2(ival)=0.0d0
          caaaa(ival)=0.0d0
          caaa2(ival)=0.0d0
          cabbb(ival)=0.0d0
          cbbbb(ival)=0.0d0
          cbbb2(ival)=0.0d0
          cbaaa(ival)=0.0d0
          caaabbb(ival)=0.0d0

          cay(ival)=0.0d0
          caaay(ival)=0.0d0
          cby(ival)=0.0d0
          cbbby(ival)=0.0d0

          inter1=interval(ival)
          inter2=interval(ival1)-1

          xfiti=xfit(ival)
          dxfit=xfit(ival1)-xfit(ival)
          h6=dxfit**2/6.0d0

          do i=inter1,inter2

            dx=xdata(i)-xfiti

            weight=1.0d0/edata(i)

            b=dx/dxfit
            a=(1.0d0-b)

            aaa=a*((a+1.0d0)*(a-1.0d0))*h6*weight
            bbb=b*((b+1.0d0)*(b-1.0d0))*h6*weight

            a=a*weight
            b=b*weight

            ca2(ival)=ca2(ival)+a*a
            cab(ival)=cab(ival)+a*b
            cb2(ival)=cb2(ival)+b*b

            caaa2(ival)=caaa2(ival)+aaa*aaa
            cbbb2(ival)=cbbb2(ival)+bbb*bbb

            caaaa(ival)=caaaa(ival)+a*aaa
            cbbbb(ival)=cbbbb(ival)+b*bbb
            cabbb(ival)=cabbb(ival)+a*bbb
            cbaaa(ival)=cbaaa(ival)+b*aaa

            caaabbb(ival)=caaabbb(ival)+aaa*bbb

            y=ydata(i)*weight

            cay(ival)=cay(ival)+a*y
            caaay(ival)=caaay(ival)+aaa*y
            cby(ival)=cby(ival)+b*y
            cbbby(ival)=cbbby(ival)+bbb*y
          enddo  !i (data)

        enddo !ival

        fmat7(1:7,nival*4)=0.0d0

        fmat7(4,1)=ca2(1)
        fmat7(5,1)=caaaa(1)
        fmat7(6,1)=cab(1)
        fmat7(7,1)=cabbb(1)
        vfit(1)=cay(1)

        fmat7(3,2)=caaaa(1)
        fmat7(4,2)=caaa2(1)
        fmat7(5,2)=cbaaa(1)
        fmat7(6,2)=caaabbb(1)
        vfit(2)=caaay(1)

        fmat7(2,3)=cab(1)
        fmat7(3,3)=cbaaa(1)
        fmat7(4,3)=cb2(1)
        fmat7(5,3)=cbbbb(1)
        vfit(3)=cby(1)

        fmat7(1,4)=cabbb(1)
        fmat7(2,4)=caaabbb(1)
        fmat7(3,4)=cbbbb(1)
        fmat7(4,4)=cbbb2(1)
        vfit(4)=cbbby(1)

        do ival=2,nival

          i4=(ival-1)*4
          i5=i4+1
          i6=i5+1
          i7=i6+1
          i8=i7+1

c last point of interval ival-1 agrees with first point of interval ival
          fmat7(1:7,i4)=0.0d0
          fmat7(3,i4)=1.0d0
          fmat7(5,i4)=-1.0d0
          vfit(i4)=0.0d0

c 2nd deriv. of last point of interval ival-1 agrees with that of first point of
c interval ival
          fmat7(1:7,i5)=0.0d0
          fmat7(5,i5)=1.0d0
          fmat7(3,i5)=-1.0d0
          vfit(i5)=0.0d0

          fmat7(3,i6)=caaaa(ival)
          fmat7(4,i6)=caaa2(ival)
          fmat7(5,i6)=cbaaa(ival)
          fmat7(6,i6)=caaabbb(ival)
          vfit(i6)=caaay(ival)

          fmat7(2,i7)=cab(ival)
          fmat7(3,i7)=cbaaa(ival)
          fmat7(4,i7)=cb2(ival)
          fmat7(5,i7)=cbbbb(ival)
          vfit(i7)=cby(ival)

          fmat7(1,i8)=cabbb(ival)
          fmat7(2,i8)=caaabbb(ival)
          fmat7(3,i8)=cbbbb(ival)
          fmat7(4,i8)=cbbb2(ival)
          vfit(i8)=cbbby(ival)

        enddo

c solve it

        call util_smooth_gauss_cs(4*nival,fmat7,vfit)

        if (ifail.ne.0) then
          cerr='*** Equation system could not be solved'
          return
        endif

c{ smoothed data

        do ival=1,nival

          ival1=ival+1
          ipar=4*(ival-1)+1

          xfiti=xfit(ival)
          dxfit=xfit(ival1)-xfit(ival)
          h6=dxfit**2/6.0d0

          do i=interval(ival),interval(ival1)-1

            dx=xdata(i)-xfiti
            b=dx/dxfit
            a=1.0d0-b
            aaa=a*(a+1.0d0)*(a-1.0d0)*h6
            bbb=b*(b+1.0d0)*(b-1.0d0)*h6

            ysmooth(i)=
     &        a*vfit(ipar)+aaa*vfit(ipar+1)+
     &        b*vfit(ipar+2)+bbb*vfit(ipar+3)

          enddo

          xcen(ival)=(xfit(ival1)+xfit(ival))/2.0d0

          a=0.5d0
          aaa=a*(a+1.0d0)*(a-1.0d0)*h6
          b=1.0d0-a
          bbb=b*(b+1.0d0)*(b-1.0d0)*h6

          ycen(ival)=
     &      a*vfit(ipar)+aaa*vfit(ipar+1)+
     &      b*vfit(ipar+2)+bbb*vfit(ipar+3)

          yfit(ival)=vfit(4*(ival-1)+1)

        enddo !ival

        yfit(nival+1)=vfit(4*(nival-1)+3)

      else if (mord.eq.4) then

        !CERN E210 package
        nw=ndatap*(nivalp+5)+nivalp*(nivalp+1)
        call dspap1(2,nknots,ndata,xdata,ydata,knotmode,xknots,coef,ws,nw,
     &    ifail)

        if (ifail.ne.0) then
          cerr='*** Error in util_smooth: Bad return from DSPAP1 ***'
          return
        endif

        if (knotmode.eq.1.or.knotmode.eq.2) then
          xfit(1:nknots)=xknots(1:nknots)
        endif

        if (kminuit.eq.0) then

          do ival=1,nival
            ivalp=ival+1
            xcen(ival)=(xfit(ivalp)+xfit(ival))/2.0d0
            ycen(ival)=dspps1(2,nknots,0,xcen(ival),xknots,coef,ifail)
            if (ifail.ne.0) then
              cerr='*** Error in util_smooth: Bad return from DSPPS1 ***'
              return
            endif
            if (ifail.ne.0) then
              cerr='*** Error in util_smooth: Bad return from DSPPS1 ***'
              return
            endif
          enddo !ival

          do ival=1,nknots
            yfit(ival)=dspps1(2,nknots,0,xfit(ival),xknots,coef,ifail)
            if (ifail.ne.0) then
              cerr='*** Error in util_smooth: Bad return from DSPPS1 ***'
              return
            endif
          enddo !nknots

        endif !(kminuit.eq.0) then

        chi2=0.0d0
        do i=1,ndata
          ysmooth(i)=dspps1(2,nknots,0,xdata(i),xknots,coef,ifail)
          if (ifail.ne.0) then
            cerr='*** Error in util_smooth: Bad return from DSPPS1 ***'
            return
          endif
          chi2=chi2+(ysmooth(i)-ydata(i))**2
        enddo

      else !mord

        cerr='*** Error in util_smooth: Bad order'
        ifail=1

        return

      endif !mord

      if (mord.ne.4) then
        chi2=0.0d0
        do ival=1,nival
          ival1=ival+1
          do i=interval(ival),interval(ival1)-1
            chi2=chi2+((ysmooth(i)-ydata(i))/edata(i))**2
          enddo
        enddo !ival
      endif

      if (kminuit.ne.0) return

      if (mdata.gt.npar) then
        chi2ndf=sqrt(chi2/(mdata-npar))
        write(6,*)
        write(6,*)' sqrt(chi2/ndf): ',chi2ndf
        write(6,*)
      else
        write(6,*)
        write(6,*)' chi2: ',chi2
        write(6,*)
      endif

      return
      end
