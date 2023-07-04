*CMZ : 00.00/15 09/10/2013  14.46.09  by  Michael Scheer
*-- Author :    Michael Scheer   09/10/2013
      subroutine util_zero(npoi,x,y,xzero,nzero)

c +PATCH,//UTIL/FOR
c +DECK,util_zero.

      implicit none

      double precision x(npoi),y(npoi),xzero(npoi)
      integer npoi,nzero,i

      nzero=0

      do i=1,npoi-1
        if (
     &      y(i).ge.0.0d0.and.y(i+1).lt.0.0d0
     &    .or.
     &      y(i).lt.0.0d0.and.y(i+1).ge.0.0d0) then
          nzero=nzero+1
          xzero(nzero)=x(i)-y(i)/(y(i+1)-y(i))*(x(i+1)-x(i))
        endif
      enddo

      return
      end
