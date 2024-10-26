*CMZ :          10/08/2024  16.46.49  by  Michael Scheer
*CMZ :  4.01/05 16/04/2024  14.41.20  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.32.18  by  Michael Scheer
*CMZ :  4.01/03 17/05/2023  10.57.05  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  13.32.32  by  Michael Scheer
*CMZ :  4.01/00 22/02/2023  14.57.49  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_fold_2d(nx,ny,x,y,fin,sigx,sigy,fold)

      implicit none

      integer ix,iy,nx,ny

      real*8 fin(nx,ny),x(nx),y(ny),sigx,sigy,fold(nx,ny),
     &  fxf(nx,ny),f(max(nx,ny)),fg(max(nx,ny)),ws1(max(nx,ny)),ws2(max(nx,ny))

      do iy=1,ny
        f(1:nx)=fin(1:nx,iy)
        if (sigx.gt.0.0d0) then
          call util_fold_function_gauss_lin(nx,x,f,sigx,3.0d0,fg,ws1,ws2)
          fxf(1:nx,iy)=fg(1:nx)
        else
          fxf(1:nx,iy)=fin(1:nx,iy)
        endif
      enddo

      do ix=1,nx
        if (sigy.gt.0.0d0) then
          f(1:ny)=fxf(ix,1:ny)
          call util_fold_function_gauss_lin(ny,y,f,sigy,3.0d0,fg,ws1,ws2)
          fold(ix,1:ny)=fg(1:ny)
        else
          fold(ix,1:ny)=fxf(ix,1:ny)
        endif
      enddo

      return
      end
