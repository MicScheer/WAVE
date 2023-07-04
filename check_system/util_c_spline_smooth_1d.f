*CMZ : 00.00/14 22/09/2011  16.08.46  by  Michael Scheer
*-- Author :    Michael Scheer   19/09/2011
      subroutine util_c_spline_smooth_1d(n,xa,ya,mknots,xknots,coef,deriv,ws,
     &  x,y,yp,ypp,yint,mode,ifail,cerr)

c simplified call to spline package NORBAS, CERN E210

c calls CERN routines in special case of interpolating cubic spline

c mode: <0 new initialization

      integer n
      double precision xa(n),ya(n),xknots(n+3),coef(n),deriv(n),
     &  ws(n*(n+5)+n*(n+1)),x,y,yp,ypp,yint,dspps1,dx

      integer ndeg,nderiv,ifail,knotmode,nknots,mode,nw,mknots,ik
      character(128) cerr

      data ndeg,nderiv,knotmode/2,1,0/

      nknots=mknots+6

      if (nknots.gt.n+ndeg+1) then
        ifail=1
        cerr='*** Error in util_c_spline_smooth_1d: Number of nodes too high ***'
        return
      endif

      if (nknots.lt.2*ndeg+2) then
        ifail=1
        cerr='*** Error in util_c_spline_smooth_1d: Number of nodes too low ***'
        return
      endif

      xknots(1:3)=xa(1)
      dx=(xa(n)-xa(1))/(n-1)
      do ik=1,mknots
        xknots(3+ik)=ik*dx
      enddo
      xknots(nknots-2:nknots)=xa(n)

      if (mode.lt.0) then
        nw=n*(nknots+5)+nknots*(nknots+1) !etwas überdimensioniert...
        call dspap1(ndeg,nknots,n,xa,ya,knotmode,xknots,coef,ws,nw,ifail)
        if (ifail.ne.0) then
          cerr='*** Error in util_c_spline_smooth_1d: Bad return from DSPIN1 ***'
          return
        endif
        call dspcd1(ndeg,nknots,nderiv,xknots,coef,deriv,ifail)
        if (ifail.ne.0) then
          cerr='*** Error in util_c_spline_smooth_1d: Bad return from DSPCD1 ***'
          return
        endif
      endif

      nderiv=0
      y    = dspps1(ndeg,nknots, nderiv,x,xknots,coef,ifail) !spline
      if (ifail.ne.0) then
        cerr='*** Error in util_c_spline_smooth_1d: Bad return from DSPPS1 for y-value ***'
        return
      endif
      nderiv=1
      yp   = dspps1(ndeg,nknots, nderiv,x,xknots,deriv,ifail) !derivative
      if (ifail.ne.0) then
        cerr='*** Error in util_c_spline_smooth_1d: Bad return from DSPPS1 for first derivative ***'
        return
      endif
      nderiv=2
      ypp   = dspps1(ndeg,nknots, nderiv,x,xknots,deriv,ifail) !derivative
      if (ifail.ne.0) then
        cerr='*** Error in util_c_spline_smooth_1d: Bad return from DSPPS1 for second derivative ***'
        return
      endif
      nderiv=-1
      yint = dspps1(ndeg,nknots, nderiv,x,xknots,coef,ifail) !integral
      if (ifail.ne.0) then
        cerr='*** Error in util_c_spline_smooth_1d: Bad return from DSPPS1 for integral ***'
        return
      endif

      return
      end
