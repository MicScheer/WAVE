*CMZ :          16/08/2025  13.56.50  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2025
      subroutine undulator_wigner_num(nx,ny,dx,dy,wlen,esour,ntx,nty,thex,they,wig,curr,banwid)

      use omp_lib

      implicit none

      include 'phyconparam.cmn'

      integer, intent(in) :: nx,ny,ntx,nty

      complex*16 :: ci=(0.0d0,1.0d0),exptx,expdtx,
     &  cthe,eki,expom,em,ep

      complex*16, intent(in) :: esour(2,nx,ny)
c      complex*16  :: wkern(2*nx,2*ny,nx,ny)

      real*8, intent(in):: wlen,thex(ntx),they(nty),dx,dy
      real*8, intent(out) :: wig(nx,ny,ntx,nty)

      real*8 :: tx,ty,ek,dtx,dty,xm,xp,xpm,ypm,ym,yp,rp,rm,x,y,wlen12,wignor,curr,banwid,
     &  specnor_si

      real secin,secout

      integer :: ix,iy,itx,ity,kx,ky,jfail,nper,iypm,ixpm,lx,ly,ifound,nmaxth=1

c      print*
c      print*,"     Calculating Wigner Distribution for ",sngl(wlen)," nm"
c      print*

c      secin=secnds(0.0)

      wlen12=1.0d0/(wlen*1.0d-9)**2
      ek=twopi1/(wlen*1.0d-9) !1/m
      eki=ci*ek

      wig=0.0d0

      if (ntx.gt.1) then
        dtx=thex(2)-thex(1)
      else
        dtx=1.0d0
      endif

      if (nty.gt.1) then
        dty=they(2)-they(1)
      else
        dty=1.0d0
      endif

      nmaxth=OMP_GET_MAX_THREADS()

c      call undulator_wigner_kernel(nx,ny,esour,wkern) bringt's nicht

      SPECNOR_SI= !merke/synchrotron_radiation.txt
     &  curr ! Strom
     &  /echarge1/hbar1*clight1/PI1*EPS01
     &  *banwid

      wignor=specnor_si*dx*dy*wlen12

!$OMP PARALLEL NUM_THREADS(nmaxth) DEFAULT(PRIVATE)
!$OMP& FIRSTPRIVATE(nx,ny,ntx,nty,dx,dy,eki,wlen12,they,thex,dtx,dty,wignor)
!$OMP& SHARED(esour,wig)

!$OMP DO

      ! Im Zentrum des Undulators, Gl. 77

      !E-Feld V/m = 1.e-4 / (clight1/1.e8) statvolt/cm) =~ 1.e-4/3.

      do iy=1,ny

        do ix=1,nx

              do iypm=-ny+1,ny-1

                ypm=dy*iypm

                ky=iy-iypm/2
                ly=iy+iypm/2

                if (ky.lt.1.or.ky.gt.ny) cycle
                if (ly.gt.ny.or.ly.lt.1) cycle

                do ity=1,nty
                  ty=they(ity)

                do ixpm=-nx+1,nx-1

                  xpm=dx*ixpm

                  kx=ix-ixpm/2
                  lx=ix+ixpm/2

                  if (kx.lt.1.or.kx.gt.nx) cycle
                  if (lx.gt.nx.or.lx.lt.1) cycle

                  do itx=1,ntx

                    em=esour(1,kx,ky)
                    ep=esour(2,lx,ly)

c                  em=wkern(ixpm+nx,iypm+ny,kx,ky)
c                  ep=wkern(ixpm+nx,iypm+ny,lx,ly)

                  if (itx.eq.1) then
                    tx=thex(itx)
                    expom=exp(-eki*(xpm*tx+ypm*ty))
                    expdtx=exp(-eki*xpm*dtx)
                  else
                    expom=expom*expdtx
                  endif

C                  expom=exp(-eki*(xpm*tx+ypm*ty))

                  wig(ix,iy,itx,ity)=wig(ix,iy,itx,ity)+
     &              dreal(em*ep*expom)*wignor

                enddo
              enddo

            enddo
          enddo

        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

c      secout=secnds(0.0)

c      print*,"     Seconds used:",secout-secin

      end
