*CMZ :  3.05/09 31/07/2018  14.23.52  by  Michael Scheer
*-- Author :    Michael Scheer   31/07/2018
      subroutine util_min_max_sngl(n,x,kmin,kmax,xmin,xmax)

      implicit none

      integer n,kmin,kmax,i
      real x(n),xmin,xmax

      xmin=1.0e30
      xmax=-1.0e30

      do i=1,n
        if(x(i).lt.xmin) then
          xmin=x(i)
          kmin=i
        endif
        if(x(i).gt.xmax) then
          xmax=x(i)
          kmax=i
        endif
      enddo

      end
