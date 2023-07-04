*CMZ : 00.00/07 29/04/2008  13.02.28  by  Michael Scheer
*CMZ : 00.00/02 25/08/2006  15.27.06  by  Michael Scheer
*CMZ : 00.00/01 23/02/96  14.56.50  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.54  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE util_index(n,xa,x,i,eps)

c returns i with:
c abs(x-xa(i)).le.eps

      implicit none

      double precision xa(n),x,eps
      integer n,i,k

      i=-1

      if (n.lt.1) return

      do k=1,n
        if (abs(x-xa(k)).lt.eps) then
          i=k
          return
        endif
      enddo

      return
      end
