*CMZ : 00.00/07 01/02/2008  16.04.42  by  Michael Scheer
*CMZ : 00.00/02 25/08/2006  15.27.06  by  Michael Scheer
*CMZ : 00.00/01 23/02/96  14.56.50  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.54  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_bisec(n,xa,x,i)

c returns i with:
c x.ge.x(i) .and. x.lt.x(i+1), if x(1).lt.x(n)
c or
c x.le.x(i) .and. x.gt.x(i+1), if x(n).lt.x(1)

      IMPLICIT NONE

      REAL*8 XA(N),x
      integer n,i,klo,khi,k

      i=-1

      if (n.lt.2) then
        return
      else if (x.eq.xa(1)) then
        i=1
        return
      else if (x.eq.xa(n)) then
        i=n-1
        return
      endif

      if (xa(1).lt.xa(n)) then

        if (x.lt.xa(1).or.x.gt.xa(n)) return

        klo=1
        KHI=N

1       IF (KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XA(K).GT.X)THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
          GOTO 1
        ENDIF

      else if (xa(n).lt.xa(1)) then

        if (x.lt.xa(n).or.x.gt.xa(1)) return

        klo=1
        Khi=n

2       IF (KHI-KLO.GT.1) THEN

          K=(KHI+KLO)/2

          IF(XA(K).LT.X)THEN
            Khi=K
          ELSE
            Klo=K
          ENDIF

          GOTO 2

        ENDIF

      endif !xa(1).lt.xa(n)

      i=klo

      RETURN
      END
