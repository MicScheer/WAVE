*CMZ :          28/10/2025  16.04.33  by  Michael Scheer
*-- Author :    Michael Scheer   28/10/2025
      subroutine util_check_grid_3d(nxyz,grid,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax)

      implicit none

      integer nxyz,nx,ny,nz,i,ix,iy,iz
      double precision grid(3,nxyz),xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,eps,x1,y1,z1,z2,dx,dy,dz
      double precision, dimension(:), allocatable :: hx,hy,hz
      integer, dimension(:), allocatable :: mx,my,mz

      nx=-1
      ny=-1
      nz=-1

      if (nxyz.lt.8) return

      xmin=1.0d30
      xmax=-1.0d30
      ymin=1.0d30
      ymax=-1.0d30
      zmin=1.0d30
      zmax=-1.0d30

      x1=grid(1,1)
      y1=grid(2,1)
      z1=grid(3,1)
      z2=grid(3,2)

      do i=1,nxyz
        x=grid(1,i)
        y=grid(2,i)
        z=grid(3,i)
        xmin=min(xmin,x)
        xmax=max(xmax,x)
        ymin=min(ymin,y)
        ymax=max(ymax,y)
        zmin=min(zmin,z)
        zmax=max(zmax,z)
      enddo

      eps=min(xmax-xmin,ymax-ymin,zmax-zmin)/1.0d15

      nz=1
      do i=2,nxyz
        z=grid(3,i)
        if (abs(z-z1).le.eps) exit
        nz=nz+1
      enddo

      ny=1
      do i=nz+1,nxyz,nz
        y=grid(2,i)
        if (abs(y-y1).le.eps) exit
        ny=ny+1
      enddo

      nx=1
      do i=nz*ny+1,nxyz,nz*ny
        x=grid(1,i)
        if (abs(x-x1).le.eps) exit
        nx=nx+1
      enddo

      if (nx*ny*nz.ne.nxyz) then
        nx=-1
        ny=-1
        nz=-1
        return
      endif

      allocate(hx(nx),hy(ny),hz(nz),mx(nx),my(ny),mz(nz))
      mx=0
      my=0
      mz=0

      dx=(xmax-xmin)/max(1,nx)
      hx(1)=xmin
      do i=2,nx
        hx(i)=hx(i-1)+dx
      enddo

      dy=(ymax-ymin)/max(1,ny)
      hy(1)=ymin
      do i=2,ny
        hy(i)=hy(i-1)+dy
      enddo

      dz=(zmax-zmin)/max(1,nz)
      hz(1)=zmin
      do i=2,nz
        hz(i)=hz(i-1)+dz
      enddo

      do i=1,nxyz
        x=grid(1,i)
        y=grid(2,i)
        z=grid(3,i)
        ix=min(int((x-xmin)/dx)+1,nx)
        mx(ix)=mx(ix)+1
        iy=min(int((y-ymin)/dy)+1,ny)
        my(iy)=my(iy)+1
        iz=min(int((z-zmin)/dz)+1,nz)
        mz(iz)=mz(iz)+1
      enddo

      do i=1,nx
        if (mx(i).ne.ny*nz) then
          nx=-1
          return
        endif
      enddo

      do i=1,ny
        if (my(i).ne.nx*nz) then
          ny=-1
          return
        endif
      enddo

      do i=1,nz
        if (mz(i).ne.nx*ny) then
          nz=-1
          return
        endif
      enddo

      deallocate(hx,hy,hz,mx,my,mz)

      return
      end
