*CMZ : 00.00/07 07/06/2011  13.46.51  by  Michael Scheer
*-- Author :    Michael Scheer   08/05/2008
      subroutine util_interpolation_parabel(n,xa,ya,x,y,yp,ypp,a3n,mode,ifail)

c *** a(3,n)!!!

      implicit none

      double precision xa(n),ya(n),x,y,yp,ypp,xa1old,xanold
      double precision a3n(3,n)

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

      klo=1

      if(mode.ge.0.and.nold.eq.n
     &    .and. xa(1).eq.xa1old
     &    .and. xa(n).eq.xanold
     &    ) then
        if (x.gt.xa(klold)) then
          klo=klold
        endif
      else
        mold=0
      endif

      if (x.lt.xa(klo+1)) then
        khi=klo+1
        goto 3
      endif

      if (xa(1).lt.xa(n)) then

        klo=1
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

        klo=1
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

      if (mold.eq.0) then

        call util_interpol_parabel_coef(n,xa,ya,a3n,ifail)
        if (ifail.ne.0) then
          return
        endif

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

      endif !mold

      y=a3n(1,klo)+a3n(2,klo)*x+a3n(3,klo)*x*x
      yp=a3n(2,klo)+2.0d0*a3n(3,klo)*x
      ypp=2.0d0*a3n(3,klo)

      ifail=0
      xa1old=xa(1)
      xanold=xa(n)
      klold=klo
      nold=n
      mold=1

      return
      end
