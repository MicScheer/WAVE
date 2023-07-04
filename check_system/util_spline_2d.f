*CMZ : 00.00/15 19/04/2013  14.50.35  by  Michael Scheer
*CMZ : 00.00/14 17/09/2011  17.12.55  by  Michael Scheer
*CMZ : 00.00/13 15/09/2011  13.20.00  by  Michael Scheer
*-- Author :    Michael Scheer   10/09/2011
      subroutine util_spline_2d(x,y,z,nx,ny,xa,ya,za,coef,coefxy,
     &  wxy,ws1,ws2,ws3,ws4,ifail,cerr)

c Ansatz: 2d the spline coefficients are splined themselve.


c +PATCH,//UTIL/FOR
c +DECK,util_spline_2d.

      implicit none

      integer matp
      parameter (matp=9) !hier

      integer nx,ny,mode,ifail,ix,iy,ical,kx,ky,
     &  kxold,kyold,klo,khi,kx1,ky1

      double precision
     &  xa(nx),ya(ny),za(nx,ny),
     &  coef(3,nx,ny),
     &  coefxy(max(nx,ny)),wxy(max(nx,ny)),
     &  ws1(max(nx,ny)),ws2(max(nx,ny)),ws3(max(nx,ny)),ws4(max(nx,ny)),
     &  x,y,z,hx,hy,f1,f2,fpp1,fpp2,hx26,hy26,ax,bx,ay,by,a31x,a31y,b31x,b31y,
     &  x1old,xnold,y1old,ynold,xold,yold

      character(128)cerr

      data ical/0/

      if (nx.le.3.or.ny.le.3) then
        ifail=1
        cerr='*** Error in util_spline_2d: nx.le.3.or.ny.le.3'
        return
      endif

      if (
     &    xa(1).ne.x1old.or.xa(nx).ne.xnold.or.
     &    ya(1).ne.y1old.or.ya(ny).ne.ynold) then
        xold=xa(nx)+1.0d30
        kxold=1
        kyold=1
        ical=0
      endif

      if (x.lt.xa(1).or.x.gt.xa(nx).or.y.lt.ya(1).or.y.gt.ya(ny)) then
        ifail=1
        cerr='*** Error in util_spline_2d: Out of range'
        z=-1.0d30
        return
      endif

      if (mode.eq.-1.or.ical.eq.0) then

        coef=-9999.0d0

        do iy=1,ny
          wxy(1:nx)=za(1:nx,iy)
          call util_spline_coef(xa,wxy,nx,9999.0d0,9999.0d0,coefxy,
     &      ws1,ws2,ws3,ws4)
          coef(1,1:nx,iy)=coefxy(1:nx)
        enddo

        do ix=1,nx
          wxy(1:ny)=za(ix,1:ny)
          call util_spline_coef(ya,wxy,ny,9999.0d0,9999.0d0,coefxy,
     &      ws1,ws2,ws3,ws4)
          coef(2,ix,1:ny)=coefxy(1:ny)
        enddo

        do ix=1,nx
          wxy(1:ny)=coef(1,ix,1:ny)
          call util_spline_coef(ya,wxy,ny,9999.0d0,9999.0d0,coefxy,
     &      ws1,ws2,ws3,ws4)
          coef(3,ix,1:ny)=coefxy(1:ny)
        enddo

      endif !mode, ical

      klo=1
      khi=nx

      if (x.ge.xa(kxold).and.x.le.xa(kxold+1)) then

        kx=kxold

      else

        if (kxold.gt.1) then
          if (x.ge.xa(kxold-1)) then
            klo=kxold-1
          endif
        endif

        if (x.ge.xa(kxold)) then
          klo=kxold
        endif

        if (kxold.lt.nx) then
          if (x.lt.xa(kxold+1)) then
            khi=kxold+1
          endif
        endif

        if (x.lt.xa(kxold)) then
          khi=kxold
        endif

1       if (khi-klo.gt.1) then
          kx=(khi+klo)/2
          if(xa(kx).gt.x)then
            khi=kx
          else
            klo=kx
          endif
          goto 1

        endif

      endif

      klo=1
      khi=ny

      if (y.ge.ya(kyold).and.y.le.ya(kyold+1)) then

        ky=kyold

      else

        if (kyold.gt.1) then
          if (y.ge.ya(kyold-1)) then
            klo=kyold-1
          endif
        endif

        if (y.ge.ya(kyold)) then
          klo=kyold
        endif

        if (kyold.gt.ny) then
          if (y.lt.ya(kyold+1)) then
            khi=kyold+1
          endif
        endif

        if (y.lt.ya(kyold)) then
          khi=kyold
        endif

11      if (khi-klo.gt.1) then
          ky=(khi+klo)/2
          if(ya(ky).gt.y)then
            khi=ky
          else
            klo=ky
          endif
          goto 11
        endif

      endif

      kx1=kx+1
      ky1=ky+1
      hx=xa(kx1)-xa(kx)
      bx=(x-xa(kx))/hx
      ax=1.0d0-bx
      hy=ya(ky1)-ya(ky)
      by=(y-ya(ky))/hy
      ay=1.0d0-by
      a31x=ax*(ax+1.0d0)*(ax-1.0d0)
      a31y=ay*(ay+1.0d0)*(ay-1.0d0)
      b31x=bx*(bx+1.0d0)*(bx-1.0d0)
      b31y=by*(by+1.0d0)*(by-1.0d0)
      hx26=hx*hx/6.0d0
      hy26=hy*hy/6.0d0

      f1=za(kx,ky)*ay+za(kx,ky1)*by+
     &  (coef(2,kx,ky)*a31y+coef(2,kx,ky1)*b31y)*hy26

      f2=za(kx1,ky)*ay+za(kx1,ky1)*by+
     &  (coef(2,kx1,ky)*a31y+coef(2,kx1,ky1)*b31y)*hy26

      fpp1=coef(1,kx,ky)*ay+coef(1,kx,ky1)*by+
     &  (coef(3,kx,ky)*a31y+coef(3,kx,ky1)*b31y)*hy26

      fpp2=coef(1,kx1,ky)*ay+coef(1,kx1,ky1)*by+
     &  (coef(3,kx1,ky)*a31y+coef(3,kx1,ky1)*b31y)*hy26

      z=f1*ax+f2*bx+(fpp1*a31x+fpp2*b31x)*hx26

      kxold=kx
      kyold=ky
      xold=x
      yold=y
      x1old=xa(1)
      xnold=xa(nx)
      y1old=ya(1)
      ynold=ya(nx)

      ical=1

      return
      end
