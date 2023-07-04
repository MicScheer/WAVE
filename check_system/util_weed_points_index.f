*CMZ :          06/06/2019  17.19.47  by  Michael Scheer
*CMZ : 00.00/20 21/01/2017  16.50.32  by  Michael Scheer
*-- Author :    Michael Scheer   20/01/2017
      subroutine util_weed_points_index(npoi,x,y,z,tolerance,iweed)

      implicit none

      double precision x(npoi),y(npoi),z(npoi),tolerance
      integer npoi,i,ifound,mpoi,k,iweed(0:npoi),iw

      mpoi=0

      do i=1,npoi

        ifound=0

        do iw=1,mpoi
          k=iweed(iw)
          if (
     &        abs(x(k)-x(i)).lt.tolerance.and.
     &        abs(y(k)-y(i)).lt.tolerance.and.
     &        abs(z(k)-z(i)).lt.tolerance) then
            ifound=k
            exit
          endif
        enddo

        if (ifound.eq.0) then
          mpoi=mpoi+1
          iweed(mpoi)=i
        endif

      enddo

      iweed(0)=mpoi

      return
      end
