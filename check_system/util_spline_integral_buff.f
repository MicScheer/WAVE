*CMZ : 00.00/16 18/03/2014  17.02.27  by  Michael Scheer
*CMZ : 00.00/07 14/08/2009  09.57.31  by  Michael Scheer
*CMZ : 00.00/02 17/08/2004  09.47.26  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_spline_integral_buff(x,y,n,resultat,nbuff,margin
     &  ,coef,work1,work2,work3,work4,istat)

c---  calculates intergral of y(x) via splines

      implicit none

      double precision x(n),y(n),resultat,
     &  coef(n),work1(n),work2(n),work3(n),work4(n),sum

      integer ibuff,n,istat,nbuff,margin,mbuff,i1,i2,nbfull

      istat=-1
      resultat=0.0d0

      mbuff=nbuff+2*margin

      if (n.le.mbuff.or.3*margin.ge.n) then
        call util_spline_integral(x,y,n,resultat
     &    ,coef,work1,work2,work3,work4)
        return
      endif

      nbfull=n/nbuff

      call util_spline_integral_window(
     &  x,y,nbuff+margin,x(1),x(nbuff),resultat
     &  ,coef,work1,work2,work3,work4,-1,istat)

      if (istat.ne.0) then
        resultat=0.0d0
        return
      endif

      i2=nbuff

      do ibuff=1,nbfull-2
        i1=i2
        i2=i1+nbuff
        call util_spline_integral_window(
     &    x(i1-margin),y(i1-margin),mbuff,x(i1),x(i2),sum
     &    ,coef,work1,work2,work3,work4,-1,istat)
        if (istat.ne.0) then
          resultat=0.0d0
          return
        endif
        resultat=resultat+sum
      enddo !nbfull

      call util_spline_integral_window(
     &  x(i2),y(i2),n-i2+1,x(i2),x(n),sum
     &  ,coef,work1,work2,work3,work4,-1,istat)

      if (istat.ne.0) then
        resultat=0.0d0
        return
      endif

      resultat=resultat+sum

      return
      end
