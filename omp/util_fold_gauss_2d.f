*CMZ :          23/08/2024  14.52.21  by  Michael Scheer
*CMZ : 00.00/15 07/12/2012  20.04.19  by  Michael Scheer
*-- Author :    Michael Scheer   06/12/2012
      subroutine util_fold_gauss_2d(nx,ny,x,y,f,sigx,rnsigx,sigy,rnsigy,fg,ispline,istat)

c Folding of f(x(ix),y(iy)) with a 2D Gaussian.
c The Gaussian is considered form -rnsig*sig -> +rnsig*sig
C IT'S ONLY CORRECT FOR X,Y FAR ENOUGH FROM THE EDGES!!

c Dimensions f(nx,ny), fg(nx,ny)

      implicit none

      double precision, dimension(:), allocatable :: wf,wfg,w1,w2,w3,w4,coef

      double precision
     &  sigx,rnsigx,sigy,rnsigy,x(nx),y(ny),f(nx,ny),fg(nx,ny)

      integer nx,ny,istat,ix,iy,ispline
      integer :: nallox=0,nalloy=0

      save

      if (2.0d0*rnsigx*sigx.ge.x(nx)-x(1).or.2.0d0*rnsigy*sigy.ge.y(ny)-y(1)) then
        istat=-1
        fg=0.0d0
        return
      endif

      if (ispline.eq.0) then
        call util_fold_gauss_lin_2d(nx,ny,x,y,f,rnsigx,sigx,rnsigy,sigy,fg)
        istat=0
        return
      endif

      if (nx.gt.0.and.istat.lt.0) deallocate(wf,wfg,w1,w2,w3,w4,coef)

      istat=0
      fg=0.0d0

      if (nx.lt.3.or.ny.lt.3) then
        istat=-1
        return
      endif

      if (nx.gt.nallox.or.ny.gt.nalloy) then
        if (nx.eq.0) deallocate(wf,wfg,w1,w2,w3,w4,coef)
        allocate(wf(max(nx,ny)))
        allocate(wfg(max(nx,ny)))
        allocate(w1(max(nx,ny)))
        allocate(w2(max(nx,ny)))
        allocate(w3(max(nx,ny)))
        allocate(w4(max(nx,ny)))
        allocate(coef(max(nx,ny)))
        nallox=nx
        nalloy=ny
      endif

      do iy=1,ny
        wf=f(1:nx,iy)
        call util_fold_function_gauss(nx,x,wf,sigx,rnsigx,wfg,coef,w1,w2,w3,w4)
        fg(1:nx,iy)=wfg(1:nx)
      enddo !iy

      do ix=1,nx
        wf=fg(ix,1:ny)
        call util_fold_function_gauss(ny,y,wf,sigy,rnsigy,wfg,coef,w1,w2,w3,w4)
        fg(ix,1:ny)=wfg(1:ny)
      enddo !iy

      return
      end
