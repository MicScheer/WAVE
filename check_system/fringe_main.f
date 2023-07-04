*CMZ :  3.05/23 12/11/2018  17.05.15  by  Michael Scheer
*-- Author :    Michael Scheer   12/11/2018
      program fringe_main

c+PATCH,//WAVE/MAIN
c +DECK,fringe_main.

      implicit none

      double precision x,y,z,bx,by,bz,ax,ay,az,fint,gap,fringe,dx,dy

      integer istatus,ix,iy,nx

      z=0.0d0
      dx=0.005d0
      nx=100
      dy=0.0010d0

      fint=0.50d0
      gap=0.025*2.0d0

      x=-dx

      do ix=1,nx
        x=x+dx
        y=-dy
        do iy=1,10
          y=y+dy
          call mrad_fringe_linear(
     &      x,y,z,bx,by,bz,ax,ay,az,fint,gap,fringe,istatus)
          write(99,*)"1. ",x,y,bx,by,istatus
        enddo
      enddo

      print*,fringe

      x=-dx

      do ix=1,nx
        x=x+dx
        y=-dy
        do iy=1,10
          y=y+dy
          call mrad_fringe_cubic_spline(
     &      x,y,z,bx,by,bz,ax,ay,az,fint,gap,fringe,istatus)
          write(99,*)"2. ",x,y,bx,by,istatus
        enddo
      enddo

      print*,fringe

      x=-dx

      do ix=1,nx
        x=x+dx
        y=-dy
        do iy=1,10
          y=y+dy
          call mrad_fringe_quintic_spline(
     &      x,y,z,bx,by,bz,ax,ay,az,fint,gap,fringe,istatus)
          write(99,*)"3. ",x,y,bx,by,istatus
        enddo
      enddo
      print*,fringe

      x=-dx

      do ix=1,nx
        x=x+dx
        y=-dy
        do iy=1,10
          y=y+dy
          call mrad_fringe_quintic_spline_msh(
     &      x,y,z,bx,by,bz,ax,ay,az,fint,gap,fringe,istatus)
          write(99,*)"4. ",x,y,bx,by,istatus
        enddo
      enddo
      print*,fringe

      end

      include 'mrad_fringe_linear.f'
      include 'mrad_fringe_cubic_spline.f'
      include 'mrad_fringe_quintic_spline.f'
      include 'mrad_fringe_quintic_spline_msh.f'
