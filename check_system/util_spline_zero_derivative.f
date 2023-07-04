*CMZ : 00.00/11 04/08/2011  15.17.40  by  Michael Scheer
*-- Author :    Michael Scheer   09/03/2011
      subroutine util_spline_zero_derivative(n,x,y,yp,ypp,xzero,yzero,ifail)

      implicit none

      integer n,ifail,null,i2,i1,i
      double precision x(n),y(n),yp(n),ypp(n),xzero,yzero,h,a,b
      double precision x1,x2,sq3,xh2,xl2,ypph2,yppl2,xh,xl,ypph,yppl,root,
     &  yl,yh,a3(3),yp3(3)

      if (n.lt.2) then
        ifail=-1
        return
      endif

      sq3=sqrt(3.0d0)

      null=0
      do i=1,n-1

        if (
     &      yp(i).le.0.0d0.and.yp(i+1).gt.0.0d0
     &      .or.
     &      yp(i).gt.0.0d0.and.yp(i+1).le.0.0d0) then

          null=null+1
          i1=i
          i2=i+1

          xl=x(i1)
          xh=x(i2)
          yl=y(i1)
          yh=y(i2)
          yppl=ypp(i1)
          ypph=ypp(i2)

          xh2=xh*xh
          xl2=xl*xl
          ypph2=ypph*ypph
          yppl2=yppl*yppl

          root=
     &      xh2*ypph2+xh2*ypph*yppl+xh2*yppl2-2.0d0*xh*xl
     &      *ypph2-2.0d0*xh*xl*ypph*yppl-2.0d0*xh*xl*yppl2+xl2*ypph2+
     &      xl2*ypph*yppl+xl2*yppl2-6.0d0*yh*ypph+6.0d0*yh*yppl+6.0d0*yl*
     &      ypph-6.0d0*yl*yppl

          if (root.gt.0.0d0) then
            root=sqrt(root)*sq3
          else
            ifail=1
            return
          endif

        endif
      enddo

      if (null.ne.1) then
        ifail=3
        return
      endif

      if (ypph-yppl.ne.0.0d0) then
        x1=(root-3.0d0*xh*yppl+3.0d0*xl*ypph)/(3.0d0*(ypph-yppl))
        x2=(-root-3.0d0*xh*yppl+3.0d0*xl*ypph)/(3.0d0*(ypph-yppl))
      else
        if (i2.lt.n) then
          call  util_parabel(x(i1),y(i1),a3,yp3,xzero,yzero,ifail)
          return
        else
          call  util_parabel(x(i2-2),y(i2-2),a3,yp3,xzero,yzero,ifail)
          return
        endif

      endif

      if (x1.ge.x(i1).and.x1.le.x(i2)) then
        xzero=x1
      else  if (x2.ge.x(i1).and.x2.le.x(i2)) then
        xzero=x2
      else
        ifail=2
        return
      endif

      xl=x(i1)
      xh=x(i2)
      yl=y(i1)
      yh=y(i2)
      h=xh-xl

      a=(xh-xzero)/h
      b=(xzero-xl)/h
      yppl=ypp(i1)
      ypph=ypp(i2)

      yzero=a*yl+b*yh+
     &  (a*(a+1.0d0)*(a-1.d0)*yppl+b*(b+1.0d0)*(b-1.0d0)*ypph)*(h**2)/6.0d0

      return
      end

