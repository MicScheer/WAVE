*CMZ :          21/01/2025  16.22.18  by  Michael Scheer
*-- Author :    Michael Scheer   20/01/2025
      module magseqf90m

      Type Magnet

      character(512) cfilefour,cfilemap,comm
      character(32) ctype,cname

      double precision, dimension (:,:), allocatable :: bmap
      double precision, dimension (:), allocatable :: four,xkfour,ykfour,zkfour

      double precision :: totlen,bmxmin,bmxmax,bmymin,bmymax,bmzmin,bmzmax,xold,yold,zold,
     &  bmbxmax,bmbxmin,bmbymax,bmbymin,bmbzmax,bmbzmin,bmapdx,bmapdy,bmapdz,a0,
     &  fouentr,fouexit,xmin,xmax,ymin,ymax,zmin,zmax,
     &  zl0four,zl0four2,xl0four,yl0four,xk0four,zk0four,yk0four,
     &  scalex=1.0d0,scaley=1.0d0,scalez=1.0d0,scalebx=1.0d0,scaleby=1.0d0,scalebz=1.0d0,
     &  offsetx=0.0d0,offsety=0.0d0,offsetz=0.0d0,offsetbx=0.0d0,offsetby=0.0d0,offsetbz=0.0d0,
     &  pmag(100),rho,rhomin,bmax,ebeam,scalebglob,xcen,ycen,zcen

      integer :: mapx=0,mapy=0,mapz=0,nfour=0,nx,ny,nz,ical=0,iwarnx=0,iwarny=0,iwarnz=0,
     &  ntot,nyz,iwarnbmap=0

      end Type Magnet

      type(Magnet), dimension(:), allocatable :: SeqMag

      end module magseqf90m
