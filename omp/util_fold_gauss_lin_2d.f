*CMZ :          16/08/2024  09.21.28  by  Michael Scheer
*-- Author :    Michael Scheer   15/08/2024
*CMZ :          15/08/2024  11.02.29  by  Michael Scheer
      subroutine util_fold_gauss_lin_2d(nx,ny,x,y,fin,rnsigx,sigx,rnsigy,sigy,fold)

      implicit none

      integer ix,iy,nx,ny

      real*8 fin(nx,ny),x(nx),y(ny),sigx,sigy,fold(nx,ny),rnsigx,rnsigy,
     &  fxf(nx,ny),f(max(nx,ny)),fg(max(nx,ny)),ws1(max(nx,ny)),ws2(max(nx,ny))

      do iy=1,ny
        f(1:nx)=fin(1:nx,iy)
        if (sigx.gt.0.0d0) then
          call util_fold_function_gauss_lin(nx,x,f,sigx,rnsigx,fg,ws1,ws2)
          fxf(1:nx,iy)=fg(1:nx)
        else
          fxf(1:nx,iy)=fin(1:nx,iy)
        endif
      enddo

      do ix=1,nx
        if (sigy.gt.0.0d0) then
          f(1:ny)=fxf(ix,1:ny)
          call util_fold_function_gauss_lin(ny,y,f,sigy,rnsigy,fg,ws1,ws2)
          fold(ix,1:ny)=fg(1:ny)
        else
          fold(ix,1:ny)=fxf(ix,1:ny)
        endif
      enddo

      return
      end
