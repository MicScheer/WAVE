*CMZ : 00.00/14 19/09/2011  15.55.03  by  Michael Scheer
*-- Author :    Michael Scheer   19/09/2011
      subroutine util_b_spline_1d(n,xa,ya,xknots,coef,deriv,ws1,iws2,
     &  x,y,yp,yint,mode,ifail,cerr)

c simplified call to spline package NORBAS, CERN E210

c calls CERN routines in special case of interpolating cubic spline

c mode: <0 new initialization

      integer n
      double precision xa(n),ya(n),xknots(n+3),coef(n),deriv(n),
     &  ws1(7*n),x,y,yp,yint,dspps1

      integer ndeg,nderiv,ifail,knotmode,nknots,mode,iws2(n)
      character(128) cerr


      data ndeg,nderiv,knotmode/2,1,1/

      nknots=n+ndeg+1

      if (mode.lt.0) then
        call dspin1(ndeg,n,xa,ya,knotmode,xknots,coef,ws1,iws2,ifail)
        if (ifail.ne.0) then
          cerr='*** Error in util_b_spline_1d: Bad return from DSPIN1 ***'
          return
        endif
        call dspcd1(ndeg,nknots,nderiv,xknots,coef,deriv,ifail)
        if (ifail.ne.0) then
          cerr='*** Error in util_b_spline_1d: Bad return from DSPCD1 ***'
          return
        endif
      endif

      nderiv=0
      y    = dspps1(ndeg,nknots, nderiv,x,xknots,coef,ifail) !spline
        if (ifail.ne.0) then
          cerr='*** Error in util_b_spline_1d: Bad return from DSPPS1 for y-value ***'
          return
        endif
      nderiv=1
      yp   = dspps1(ndeg,nknots, nderiv,x,xknots,deriv,ifail) !derivative
        if (ifail.ne.0) then
          cerr='*** Error in util_b_spline_1d: Bad return from DSPPS1 for first derivative ***'
          return
        endif
      nderiv=-1
      yint = dspps1(ndeg,nknots, nderiv,x,xknots,coef,ifail) !integral
        if (ifail.ne.0) then
          cerr='*** Error in util_b_spline_1d: Bad return from DSPPS1 for integral ***'
          return
        endif

      return
      end
