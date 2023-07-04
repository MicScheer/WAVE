*CMZ :  3.05/22 08/11/2018  13.48.34  by  Michael Scheer
*-- Author :    Michael Scheer   08/11/2018
      subroutine util_interpol_extrapol_on_grid(
     &  ndat,xdat,ydat,ngrid,xgrid,ygrid,istat)

c +PATCH,//UTIL/MAIN
c +DECK,util_interpol_extrapol_on_grid.

      implicit none

      double precision xdat(ndat),ydat(ndat),xgrid(ngrid),ygrid(ngrid)

      integer ndat,ngrid,i,istat,ifail

      do i=1,ngrid
        call util_interpol_extrapol_linear(ndat,xdat,ydat,xgrid(i),ygrid(i),ifail)
        istat=istat+ifail
      enddo

      return
      end
