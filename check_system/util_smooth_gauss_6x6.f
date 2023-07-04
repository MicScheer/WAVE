*CMZ : 00.00/11 17/08/2011  11.41.57  by  Michael Scheer
*-- Author :    Michael Scheer   17/08/2011
      subroutine util_smooth_gauss_6x6(n,a,b)

c +PATCH,//UTIL/FOR
c +DECK,util_smooth_gauss_6x6.

c Solves linear equation system for a second order spline fit matrix.
c To reduce storage size by using util_smooth_gauss_6x6_s

      integer n,i1,i2,ne

      double precision a(n,n),b(n),eps

      data eps/1.0d-20/

      do i1=1,n

        ne=min(n,i1+3)
        b(i1)=b(i1)/a(i1,i1)
        a(i1,i1:ne)=a(i1,i1:ne)/a(i1,i1)

        do i2=i1+1,ne
          if (abs(a(i2,i1)).ge.eps) then
            b(i2)=b(i2)/a(i2,i1)-b(i1)
            a(i2,i1:ne)=a(i2,i1:ne)/a(i2,i1)-a(i1,i1:ne)
          endif
        enddo

      enddo

c now we have a triangle matrix with
c      a(1,1)=1, a(2:3).ne.0
c      a(i,i)=1. and a(2:n-1,i+1).ne.0 for i=2:n-1
c      a(n:n)=1.0

      do i1=n,3,-1
        i2=i1-1
        b(i2)=(b(i2)/a(i2,i1)-b(i1))/(a(i2,i2)/a(i2,i1))
      enddo

      b(1)=(b(1)-a(1,2)*b(2)-a(1,3)*b(3))/a(1,1)

      return
      end
