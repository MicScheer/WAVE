*CMZ :          09/01/2019  13.49.55  by  Michael Scheer
*CMZ : 00.00/20 13/10/2016  12.39.38  by  Michael Scheer
*-- Author :    Michael Scheer   11/10/2016
      subroutine util_fit_2d_parabola(npoi,x,y,f,param,xopt,yopt,fopt,ifail)

c +PATCH,//UTIL/FOR
c +DECK,util_fit_2d_parabola.

      implicit none

      integer ifail,npar,npoi,narg,nfun,ndimpoi,ipoi
      parameter (npar=6,nfun=1,narg=2,ndimpoi=1000)

      real*8 fundata(narg+nfun,ndimpoi),x(npoi),y(npoi),f(npoi),det
      real*8 param(npar),a(npar,npar),t(nfun,npar)
      double precision xopt,yopt,fopt

      if (npoi.gt.ndimpoi) then
        ifail=1
        return
      endif

      do ipoi=1,npoi
        fundata(1,ipoi)=x(ipoi)
        fundata(2,ipoi)=y(ipoi)
        fundata(3,ipoi)=f(ipoi)
      enddo

      call util_linear_fit
     &  (ifail,npar,param,ndimpoi,npoi,narg,nfun,a,t,fundata)

      if (ifail.ne.0) then
        ifail=2
        return
      endif

c search extremum (x0,y0)

c 2*x0*param(4)+  y0*param(6) = -param(2)
c   x0*param(6)+2*y0*param(5) = -param(3)

      det=param(4)*param(5)-param(6)**2

      if (det.ne.0.0d0) then
        xopt=(param(2)*param(5)-param(3)*param(6))/det
        yopt=(param(4)*param(3)-param(6)*param(2))/det
        fopt=param(1)+param(2)*xopt+param(3)*yopt+param(4)*xopt*xopt+param(5)*yopt*yopt+param(6)*xopt*yopt
      else
        ifail=3
      endif

      return
      end

      include 'util_linear_fit.f'
      include 'util_fit_2d_parabola_user.f'
