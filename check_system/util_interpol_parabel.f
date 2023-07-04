*CMZ : 00.00/07 09/05/2008  08.29.47  by  Michael Scheer
*-- Author :    Michael Scheer   08/05/2008
      subroutine util_interpol_parabel(n,xa,ya,x,y,yp,ypp,mode,ifail)

c Direct mode, i.e. parabola is calculated for current interval only
c See also util_interpolation_parabel and util_interpol_parabel_coef

      implicit none

      double precision xa(n),ya(n),x,y,yp,ypp,xa1old,xanold
      double precision a(3),dxm,dxp,x0,a1,a2,dxm2,dxp2
      double precision det,a22,fm,fp,f0,dx

      integer nold,n,klo,khi,klold,k,mode,klm,klp,mold,ifail

      data klold/1/,nold/-99/
      data mold/0/
      data xa1old/-9999.0d0/,xanold/-9999.0d0/

      ifail=-1
      if (mode.lt.0) mold=0

      if (n.lt.1) then
        mold=0
        return
      else if (n.lt.3) then
        mold=0
        if (xa(1).eq.xa(2)) then
          return
        else
          ypp=0.0d0
          yp=(ya(2)-ya(1))/(xa(2)-xa(1))
          y=ya(1)+yp*(x-xa(1))
        endif
        return
      endif

      if (mode.lt.0.or.klold.ge.n) then
        mold=0
        klo=1
      else if(nold.eq.n
     &    .and. xa(1).eq.xa1old
     &    .and. xa(n).eq.xanold
     &    .and. x.gt.xa(klold)
     &    ) then
        klo=klold
      else
        mold=0
        klo=1
      endif

      if (x.lt.xa(klo+1)) then
        khi=klo+1
        goto 3
      endif

      if (xa(1).lt.xa(n)) then

        khi=n

1       if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if(xa(k).gt.x)then
            khi=k
          else
            klo=k
          endif
          goto 1
        endif

      else if (xa(n).lt.xa(1)) then

        if (x.lt.xa(n).or.x.gt.xa(1)) return

        khi=n

2       if (khi-klo.gt.1) then

          k=(khi+klo)/2

          if(xa(k).lt.x)then
            khi=k
          else
            klo=k
          endif

          goto 2

        endif

      endif !xa(1).lt.xa(n)

3     continue

      if (mold.eq.0.or.klo.ne.klold) then

        if (klo.gt.1.and.klo.lt.n) then
          klm=klo-1
          klp=klo+1
        else if (klo.eq.1) then
          klm=klo
          klo=klm+1
          klp=klo+1
        else if (klo.eq.n) then
          klp=n
          klo=klp-1
          klm=klo-1
        endif !klo

c calculate f=a0+a1*(x-x0)+a2*(x-x0)**2
c  = a0 + a1*x - a0*x0 + a2*x**2 - 2*a2*x*x0 + a2*x0**2
c  = a0 + (a2*x0 -a0)*x0 + (a1 - 2*a2*x0 )*x+ a2*x**2

c change system: (x0,s0)->(0,0), i.e.
c calculate f=a1*dx+a2*dx**2
c  df/dx=a1+2*a2*dx_max =! 0, dx_max=-a1/2/a2

        x0=xa(klo)
        f0=ya(klo)

        fm=ya(klm)-f0
        fp=ya(klp)-f0

        dxm=xa(klm)-x0
        dxp=xa(klp)-x0

c fm=a1*dxm+a2*dxm**2
c fp=a1*dxp+a2*dxp**2

c (dxm dxm2) (a1) = (y(1))
c (dxp dxp2) (a2) = (y(3))

        dxm2=dxm*dxm
        dxp2=dxp*dxp

        det=dxm*dxp2-dxp*dxm2

        if (det.ne.0.0d0) then

          a1=(fm*dxp2-fp*dxm2)/det
          a2=(fp*dxm-fm*dxp)/det

c calculate f=      f0      +  a1*dx  + a2*dx**2
c            = a1*x - a1*x0 + a2*x**2 + a2*x0**2 - 2*a2*x*x0
c  f = f0 + (a2*x0 -a1)*x0 + (a1 - 2*a2*x0 )*x+ a2*x**2

          a22=2.0d0*a2

          a(1)=f0 + (a2*x0 -a1)*x0
          a(2)=a1 - a22*x0
          a(3)=a2

c calculate yp=a1+2*a2*dx

        endif !det

      endif !mold

      dx=x-x0

      y=f0+a1*dx+a2*dx*dx
      yp=a1+2.0d0*a2*dx
      ypp=2.0d0*a2

      ifail=0
      nold=n
      klold=klo
      xa1old=xa(1)
      xaNold=xa(n)
      mold=1

      return
      end
