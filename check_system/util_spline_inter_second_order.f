*CMZ : 00.00/11 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/09 07/06/2011  13.46.51  by  Michael Scheer
*CMZ : 00.00/07 08/09/2009  09.23.13  by  Michael Scheer
*CMZ : 00.00/05 06/03/2007  16.31.51  by  Michael Scheer
*CMZ : 00.00/02 25/08/2006  15.27.06  by  Michael Scheer
*CMZ : 00.00/01 23/02/96  14.56.50  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.54  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_spline_inter_second_order(xa,ya,yap,n,x,y,yp,ypp,mode,istat)

c---  interpolates y(x) via second order spline

c--   input:

c-       n: number of x,y-values
c-       xa:   array of x-values
c-       ya:   array of y-values
c-       yap:  array spline coefficients, i.e. first derivatives
c-       x: y(x), yp(x), and yp(x) are calculated
c-       mode: control flag:
c-             mode.ge.0: use values of last call to start with
c-             mode.lt.0: new initialization

c--   output:

c-       y: y(x),yp(x), and ypp(x) are calculated

      implicit none

      integer nold,n,klo,khi,klold,k,mode,istat

      real*8 y,x,xa1old,xanold,h,eps,yp

      real*8 xa(n),ya(n),yap(n),xx,ypp

      data klold/1/,nold/-99/
      data xa1old/-9999.d0/,xanold/-9999./

      istat=0

      eps=abs(xa(n)-xa(1))/1.0d10
      xx=x

      if(xa(1).lt.xa(n)) then

        if(xx.lt.xa(1).and.xx.gt.xa(1)-eps) then
          xx=xa(1)
        else if(xx.gt.xa(n).and.xx.lt.xa(n)+eps) then
          xx=xa(n)
        endif

        if(xx.lt.xa(1).or.xx.gt.xa(n)) then
          write(6,*)'xa(1), xa(n):',xa(1), xa(n)
          write(6,*)'x:',x
          write(6 ,*)'***error in util_spline_inter: x out of range ***'
          stop
        endif

      else

        if(xx.lt.xa(n).and.xx.gt.xa(n)-eps) then
          xx=xa(n)
        else if(xx.gt.xa(1).and.xx.lt.xa(n)+eps) then
          xx=xa(1)
        endif

        if(xx.lt.xa(n).or.xx.gt.xa(1)) then
          write(6,*)'xa(1), xa(n):',xa(1), xa(n)
          write(6,*)'x:',x
          write(6 ,*)'***error in util_spline_inter: x out of range ***'
          istat=-1
          return
        endif

      endif

      if (mode.lt.0.or.klold.ge.n) then
        klo=1
      else if(nold.eq.n
     &    .and. xa(1).eq.xa1old
     &    .and. xa(n).eq.xanold
     &    .and. xx.gt.xa(klold)
     &    ) then
        klo=klold
      else
        klo=1
      endif

      if (xx.lt.xa(klo+1)) then
        khi=klo+1
        goto 2
      endif

      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.xx)then
          khi=k
        else
          klo=k
        endif
        goto 1
      endif

2     h=xa(khi)-xa(klo)

      if (h.le.0.0d0) then
        istat=-2
        return
      endif

      ypp=2.0d0*(ya(khi)-ya(klo)-yap(klo)*h)/h**2
      yp=yap(klo)+ypp*(xx-xa(klo))
      y=ya(klo)+yap(klo)*(xx-xa(klo))+ypp/2.0d0*(xx-xa(klo))**2

      klold=klo
      nold=n
      xa1old=xa(1)
      xanold=xa(n)

      return
      end
