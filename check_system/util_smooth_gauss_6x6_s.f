*CMZ : 00.00/11 19/08/2011  12.09.31  by  Michael Scheer
*-- Author :    Michael Scheer   17/08/2011
      subroutine util_smooth_gauss_6x6_s(n,a,b)

c +PATCH,//UTIL/FOR
c +DECK,util_smooth_gauss_6x6_s.

c Solves linear equation system for a second order spline fit matrix.
c Uses 6 x n array instead of n x n array as in util_smooth_gauss_6x6.f
c The index a(i,i) of the virtual n x n matrix corresponds to a(3,i) of the
c 6 x n matrix

      integer n,i1,i2

      double precision a(6,n),b(n),eps,as
      data eps/1.0d-20/

      do i1=1,n

        as=a(3,i1)
        b(i1)=b(i1)/as
        a(3:5,i1)=a(3:5,i1)/as

        i2=i1+1
        if (i2.le.n) then
          as=a(2,i2)
          if (abs(as).ge.eps) then
            b(i2)=b(i2)/as-b(i1)
            a(2:5,i2)=a(2:5,i2)/as-a(3:6,i1)
          endif
        endif

        i2=i1+2
        if (i2.le.n) then
          as=a(1,i2)
          if (abs(as).ge.eps) then
            b(i2)=b(i2)/as-b(i1)
            a(1:4,i2)=a(1:4,i2)/as-a(3:6,i1)
          endif
        endif


      enddo

      do i1=n,3,-1
        i2=i1-1
        b(i2)=(b(i2)/a(4,i2)-b(i1))/(a(3,i2)/a(4,i2))
      enddo

      b(1)=b(1)-a(4,1)*b(2)-a(5,1)*b(3)

      return
      end
