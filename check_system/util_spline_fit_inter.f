*CMZ : 00.00/19 11/08/2015  16.20.16  by  Michael Scheer
*CMZ : 00.00/18 11/08/2015  15.40.57  by  Michael Scheer
*CMZ : 00.00/17 26/07/2015  15.17.12  by  Michael Scheer
*CMZ : 00.00/16 24/07/2015  13.16.38  by  Michael Scheer
      subroutine util_spline_fit_inter(n,xa,fitpar,x,y,yp,ypp,ifail,mode)

c +PATCH,//UTIL/FOR
c +DECK,util_spline_fit_inter.

      implicit none

      integer n,ifail,mode,klo,khi,klold,nold,ipar,k
      double precision xa(n),fitpar(2*n),x,y,yp,ypp,xold,xa1old,xanold,
     &  h26,h,a3a,b3b,a2,b2,a,b

      save

      data xold/-1.0d30/
      data klold/-1/
      data xa1old/-1.0d30/
      data xanold/1.0d30/

      ifail=-1

      if (mode.lt.0. or. klold.ge.n .or. klold.lt.0) then
        klo=1
      else if(nold.eq.n
     &    .and. xa(1).eq.xa1old
     &    .and. xa(n).eq.xanold
     &    .and. x.gt.xa(klold)
     &    ) then
        klo=klold
      else
        klo=1
      endif

      if (x.lt.xa(klo+1)) then
        khi=klo+1
        goto 2
      endif

      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
        goto 1
      endif

2     h=xa(khi)-xa(klo)

      if (h.le.0.0d0) return

      h26=h*h/6.0d0
      a=(xa(klo+1)-x)/h
      b=1.0d0-a
      a2=a*a
      b2=b*b
      a3a=a*(a+1.0d0)*(a-1.0d0)*h26
      b3b=b*(b+1.0d0)*(b-1.0d0)*h26

      y=0.0d0
      yp=0.0d0

      ipar=1+(klo-1)*2

      y=y+fitpar(ipar)*a
      yp=yp-fitpar(ipar)/h

      ipar=ipar+1
      y=y+fitpar(ipar)*a3a
      yp=yp-fitpar(ipar)*(3.d0*a2-1.0d0)/6.0d0*h

      ipar=ipar+1
      y=y+fitpar(ipar)*b
      yp=yp+fitpar(ipar)/h

      ipar=ipar+1
      y=y+fitpar(ipar)*b3b
      yp=yp+fitpar(ipar)*(3.d0*b2-1.0d0)/6.0d0*h

      ipar=(klo-1)*2+2
      ypp=fitpar(ipar)+(fitpar(ipar+2)-fitpar(ipar))*b

      nold=n
      klold=klo
      xold=x
      ifail=0

      return
      end
