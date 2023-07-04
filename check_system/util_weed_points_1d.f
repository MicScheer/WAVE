*CMZ :          30/05/2019  21.07.37  by  Michael Scheer
*CMZ : 00.00/20 07/02/2017  14.27.25  by  Michael Scheer
*CMZ :  1.11/04 21/01/2017  16.50.32  by  Michael Scheer
*-- Author :    Michael Scheer   20/01/2017
      subroutine util_weed_points_1d(npoi,x,tolerance)

      implicit none

      double precision x(npoi),tolerance
      integer npoi,i,ifound,mpoi,k

      mpoi=0
      do i=1,npoi
        ifound=0
        do k=1,mpoi
          if (abs(x(k)-x(i)).lt.tolerance) then
            ifound=k
            exit
          endif
        enddo
        if (ifound.eq.0) then
          mpoi=mpoi+1
          x(mpoi)=x(i)
        endif
      enddo

      npoi=mpoi

      return
      end
