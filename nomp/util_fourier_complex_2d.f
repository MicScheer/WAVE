*CMZ :  4.01/05 28/03/2024  18.50.01  by  Michael Scheer
*CMZ : 00.00/02 22/04/99  17.27.52  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_fourier_complex_2d(nx,ny,x,y,f,nomx,nomy,omx,omy,ft,ispline,ifail)

      implicit none

      integer nx,ny,nomx,nomy,ispline,ifail

      complex*16 :: f(nx,ny),ft(nomx,nomy)

      real*8 x(nx),y(ny),omx(nomx),omy(nomy)

      if (ispline.eq.0) then
        call util_fourier_spline_complex_2d(nx,ny,x,y,f,nomx,nomy,omx,omy,ft,ifail)
      else
        call util_fourier_linear_complex_2d(nx,ny,x,y,f,nomx,nomy,omx,omy,ft,ifail)
      endif

      return
      end
