*CMZ :  3.05/22 08/11/2018  13.08.23  by  Michael Scheer
*-- Author :    Michael Scheer   08/11/2018
      subroutine util_interpol_on_grid(ndat,xdat,ydat,ngrid,xgrid,ygrid,istat)

c +PATCH,//UTIL/MAIN
c +DECK,util_interpol_on_grid.

      implicit none

      double precision xdat(ndat),ydat(ndat),xgrid(ngrid),ygrid(ngrid)

      integer ndat,ngrid,i,ifail,istat

      do i=1,ngrid
        call util_interpol_linear(ndat,xdat,ydat,xgrid(i),ygrid(i),ifail)
        if (ifail.ne.0) then
          print*,"*** Error, Bad return from util_interpol_linear ***"
          istat=i
          return
        endif
      enddo

      return
      end
