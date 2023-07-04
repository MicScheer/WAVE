*CMZ : 00.00/07 29/04/2008  12.48.09  by  Michael Scheer
*CMZ : 00.00/02 25/08/2006  15.27.06  by  Michael Scheer
*CMZ : 00.00/01 23/02/96  14.56.50  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.54  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE util_index_bisec(n,xa,x,i,eps)

c returns i with:
c abs(x-xa(i)).le.eps

c xa must be monoton

      implicit none

      double precision xa(n),x,eps
      integer n,i,klo,khi,k

      i=-1

      if (n.lt.1) then
        return
      else if (n.eq.1) then
        if (abs(x-xa(1)).lt.eps) i=1
        return
      else if (abs(x-xa(n)).lt.eps) then
        i=n
        return
      else if (abs(x-xa(n-1)).lt.eps) then
        i=n-1
        return
      endif

      if (xa(1).lt.xa(n)) then

        if (x.lt.xa(1).or.x.gt.xa(n)) return

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

      if (abs(x-xa(klo)).lt.eps) then
        i=klo
      else if (abs(x-xa(khi)).lt.eps) then
        i=khi
      endif

      return
      end
