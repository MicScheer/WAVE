*CMZ : 00.00/16 10/10/2014  16.14.26  by  Michael Scheer
*CMZ : 00.00/14 16/09/2011  15.02.41  by  Michael Scheer
*CMZ : 00.00/12 01/09/2011  11.09.00  by  Michael Scheer
      subroutine util_smooth_2d_gauss(n,a,b)

c +PATCH,//UTIL/FOR
c +DECK,util_smooth_2d_gauss.

c Solves linear equation system for a 2d cubic spline fit matrix.
c The index a(i,i) of the virtual n x n matrix corresponds to a(16,i) of the
c 32 x n matrix

      integer n,i1,i2,ne,i12,i162

      double precision a(32,n),b(n),eps,as
      data eps/1.0d-20/

      do i1=1,n
        as=a(16,i1)
        if (abs(as).gt.eps) then
          b(i1)=b(i1)/as
          a(1:32,i1)=a(1:32,i1)/as
          ne=min(i1+31,n)
          do i2=i1+1,ne
            i12=max(1,16-i2+i1)
            as=a(i12,i2)
            if (abs(as).gt.eps) then
              b(i2)=b(i2)/as
              a(i12:i12+15,i2)=a(i12:i12+15,i2)/as
              b(i2)=b(i2)-b(i1)
              a(i12:i12+15,i2)=a(i12:i12+15,i2)-a(17:32,i1)
            endif
          enddo
        endif
      enddo !n

      i1=n
      as=a(16,i1)
      if (abs(as).gt.eps) then
        b(i1)=b(i1)/as
        a(1:32,i1)=a(1:32,i1)/as
      endif

c now the matrix is a triangle, last row already done, i.e. a(16,n)=1

      do i1=n,1,-1
        b(i1)=b(i1)/a(16,i1)
        ne=min(15,i1-1)
        do i2=1,ne
          i12=i1-i2
          i162=16+i2
          as=a(i162,i12)
          if (abs(as).gt.eps) then
            a(16:i162-1,i12)=a(16:i162-1,i12)/as
            b(i12)=b(i12)/as-b(i1)
          endif
        enddo
      enddo !n

      return
      end
