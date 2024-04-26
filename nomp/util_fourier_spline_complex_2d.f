*CMZ :  4.01/05 28/03/2024  18.45.57  by  Michael Scheer
*CMZ : 00.00/02 22/04/99  17.27.52  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_fourier_spline_complex_2d(nx,ny,x,y,f,nomx,nomy,omx,omy,ft,ifail)

      implicit none

      integer nx,ny,nomx,nomy,ix,iy,ifail,jfail

      complex*16 :: f(nx,ny),ft(nomx,nomy),fxy(max(nx,ny)),ftxy(max(nx,ny)),ftx(nomx,ny)

      real*8 x(nx),y(ny),omx(nomx),omy(nomy)

      ifail=0

      do iy=1,ny
        fxy(1:nx)=f(1:nx,iy)
        call util_fourier_spline_complex(nx,x,fxy,nomx,omx,ftxy,jfail)
        ftx(1:nomx,iy)=ftxy(1:nomx)
        ifail=ifail+jfail
      enddo !ny

      do ix=1,nomx
        fxy(1:ny)=ftx(ix,1:ny)
        call util_fourier_spline_complex(ny,y,fxy,nomy,omy,ftxy,jfail)
        ft(ix,1:nomy)=ftxy(1:nomy)
        ifail=ifail+jfail
      enddo !ny

      return
      end
