*CMZ :  3.05/23 23/11/2018  18.06.23  by  Michael Scheer
*CMZ :  1.16/00 07/05/2017  12.11.52  by  Michael Scheer
*CMZ :  1.15/10 10/04/2017  12.35.04  by  Michael Scheer
*CMZ : 00.00/15 15/04/2013  19.27.04  by  Michael Scheer
*CMZ : 00.00/07 09/05/2008  08.29.47  by  Michael Scheer
*-- Author :    Michael Scheer   08/05/2008
      subroutine util_interpol_extrapol_linear(n,xa,ya,x,y,ifail)

      implicit none

      double precision xa(n),ya(n),x,y,xa1old,xanold

      integer nold,n,klo,khi,klold,k,ifail

      save

      data klold/1/,nold/-99/
      data xa1old/-9999.0d0/,xanold/-9999.0d0/

      ifail=-1
      if (n.lt.2) return

      if (xa(1).eq.xa(2).or.xa(n-1).eq.xa(n)) then
        return
      else if (xa(1).lt.xa(2)) then
        if (x.lt.xa(1)) then
          y=ya(1)+(ya(2)-ya(1))/(xa(2)-xa(1))*(x-xa(1))
          ifail=0
          return
        else if (x.gt.xa(n)) then
          y=ya(n)+(ya(n)-ya(n-1))/(xa(n)-xa(n-1))*(x-xa(n))
          ifail=0
          return
        endif
      else
        if (x.gt.xa(1)) then
          y=ya(1)+(ya(2)-ya(1))/(xa(2)-xa(1))*(x-xa(1))
          ifail=0
          return
        else if (x.lt.xa(n)) then
          y=ya(n)+(ya(n)-ya(n-1))/(xa(n)-xa(n-1))*(x-xa(n))
          ifail=0
          return
        endif
      endif

      klo=1
      khi=n

      klo=1

      if (klold.ge.n) then
        khi=n
      else if(nold.eq.n
     &    .and. xa(1).eq.xa1old
     &    .and. xa(n).eq.xanold
     &    ) then
        if (x.gt.xa(klold)) then
          klo=klold
        else if (klold.gt.1) then
          if (x.gt.xa(klold-1)) klo=klold-1
        endif
      endif

      if (x.lt.xa(klo+1)) then
        khi=klo+1
        goto 3
      endif

      if (xa(1).lt.xa(n)) then

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

3     continue

      y=ya(klo)+(ya(khi)-ya(klo))/(xa(khi)-xa(klo))*(x-xa(klo))

      ifail=0
      nold=n
      klold=klo
      xa1old=xa(1)
      xanold=xa(n)

      return
      end
