*CMZ : 00.00/15 26/10/2012  09.35.15  by  Michael Scheer
*-- Author :    Michael Scheer   25/10/2012
      subroutine util_straight_line_fit(n,x,y,e,a,b,chi2,erra,errb,istat)

c Fit a and b, such that f = a*x+b and Sum(i,((a*xi+b-yi-fi)/ei)**2)=min
c See also grafit.kumac

c The uncertainties erra and errb of the fitted parameters a (slope) and b
c (offset) agree to the fits of MINUIT, expept when the error e(i) are zero.
c Then the errors are such, that chi2/ndf is scaled to unity

      integer istat,n,i

      double precision x(n),y(n),e(n),chi2,wsum,wxy,wx,wy,wx2,
     &  d,ep,erra,errb,a,b,w,esuma

      wx=0.0d0
      wy=0.0d0
      wxy=0.0d0
      wx2=0.0d0
      wsum=0.0d0
      chi2=0.0d0
      erra=0.0d0
      errb=0.0d0
      istat=0

      if (n.lt.2) then
        istat=-1
        return
      else if (n.eq.2) then
        if (x(2).eq.x(1)) then
          istat=-2
          return
        endif
        a=(y(2)-y(1))/(x(2)-x(1))
        b=y(2)-a*x(2)
        return
      endif

      ep=1.0d0
      esuma=0.0d0
      do i=1,n
        if (e(i).lt.0.0d0) then
          istat=1
          return
        else if (e(i).eq.0.0d0) then
          ep=0.0d0
        endif
        esuma=esuma+abs(e(i))
      enddo

      if(esuma.ne.0.0d0.and.ep.eq.0.0d0) then
        istat=2
        return
      endif

      do i=1,n
        if (esuma.ne.0.0d0) then
          w=1.0d0/e(i)**2
        else
          w=1.0d0
        endif
        wx=wx+w*x(i)
        wy=wy+w*y(i)
        wx2=wx2+w*x(i)*x(i)
        wxy=wxy+w*x(i)*y(i)
        wsum=wsum+w
      enddo

      d=wsum*wx2-wx**2

      if (d.ne.0.0d0) then

        a=(wsum*wxy-wx*wy)/d
        b=(wx2*wy-wx*wxy)/d

        if (esuma.ne.0.0d0) then
          do i=1,n
            chi2=chi2+((a*x(i)+b-y(i))/e(i))**2
          enddo
          erra=sqrt(wsum/d)
          errb=sqrt(wx2/d)
        else
          do i=1,n
            chi2=chi2+(a*x(i)+b-y(i))**2
          enddo
          erra=sqrt(wsum/d*chi2/(n-2))
          errb=sqrt(wx2/d*chi2/(n-2))
        endif

      else
        istat=-1
        return
      endif

      return
      end
