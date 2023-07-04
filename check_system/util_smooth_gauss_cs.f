*CMZ : 00.00/12 01/09/2011  11.09.00  by  Michael Scheer
*CMZ : 00.00/11 19/08/2011  12.09.31  by  Michael Scheer
      subroutine util_smooth_gauss_cs(n,a,b)

c +PATCH,//UTIL/FOR
c +DECK,util_smooth_gauss_cs.

c Solves linear equation system for a cubic spline fit matrix.
c The index a(i,i) of the virtual n x n matrix corresponds to a(4,i) of the
c 7 x n matrix


c x x x x | 0 0
c x x x x | 0 0
c x x x x | 0 0
c 1 0 0 0 | 1 0
c--------------
c 0 1 0 0 | 0 1

      integer n,i1,i2,ne,i12,i42

      double precision a(7,n),b(n),eps,as
      data eps/1.0d-20/

      do i1=1,n-1
        as=a(4,i1)
        if (abs(as).gt.eps) then
          b(i1)=b(i1)/as
          a(1:7,i1)=a(1:7,i1)/as
          ne=min(i1+3,n)
          do i2=i1+1,ne
            i12=max(1,4-i2+i1)
            as=a(i12,i2)
            if (abs(as).gt.eps) then
              b(i2)=b(i2)/as
              a(i12:i12+3,i2)=a(i12:i12+3,i2)/as
              b(i2)=b(i2)-b(i1)
              a(i12:i12+3,i2)=a(i12:i12+3,i2)-a(4:7,i1)
            endif
          enddo
        endif
      enddo !n

      i1=n
      as=a(4,i1)
      if (abs(as).gt.eps) then
        b(i1)=b(i1)/as
        a(1:7,i1)=a(1:7,i1)/as
      endif
c now the matrix is a triangle, last row already done, i.e. a(4,n)=1

      do i1=n,1,-1
        b(i1)=b(i1)/a(4,i1)
        ne=min(3,i1-1)
        do i2=1,ne
          i12=i1-i2
          i42=4+i2
          as=a(i42,i12)
          if (abs(as).gt.eps) then
            a(4:i42-1,i12)=a(4:i42-1,i12)/as
            b(i12)=b(i12)/as-b(i1)
          endif
        enddo
      enddo !n

      return
      end
