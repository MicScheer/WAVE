*CMZ :          04/10/2025  07.34.01  by  Michael Scheer
*CMZ :  4.02/00 11/09/2025  15.05.09  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2025
      subroutine undulator_wigner_num(nx,ny,dx,dy,wlen,esour,ntx,nty,thex,they,wig,curr,banwid)

      use omp_lib

      implicit none

      include 'phyconparam.cmn'

      integer, intent(in) :: nx,ny,ntx,nty

      complex*16 :: ci=(0.0d0,1.0d0),exptx,expdtx,
     &  cthe,eki,expom,em,ep

      complex*16, intent(in) :: esour(2,nx,ny)

      real*8, intent(in):: wlen,thex(ntx),they(nty),dx,dy
      real*8, intent(out) :: wig(nx,ny,ntx,nty)

      real*8 :: tx,ty,ek,dtx,dty,xm,xp,xpm,ypm,ym,yp,rp,rm,x,y,wlen12,wignor,curr,banwid,
     &  specnor_si

      integer :: ix,iy,itx,ity,kx,ky,nper,iypm,ixpm,lx,ly,nmaxth=1

      wlen12=1.0d0/(wlen*1.0d-9)**2
      ek=twopi1/abs(wlen*1.0d-9) !1/m
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

      SPECNOR_SI= !merke/synchrotron_radiation.txt
     &  curr ! Strom
     &  /echarge1/hbar1*clight1/PI1*EPS01
     &  *banwid

      wignor=specnor_si*dx*dy*wlen12 !??*4.0d0 ! factor 4 due to change in integration variables!?

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

                  if (wlen.gt.0.0d0) then
                    if (itx.eq.1) then
                      tx=thex(itx)
                      expom=exp(-eki*(xpm*tx+ypm*ty))
                      expdtx=exp(-eki*xpm*dtx)
                    else
                      expom=expom*expdtx
                    endif
                  else
                    if (itx.eq.1) then
                      tx=thex(itx)
                      expom=exp(eki*(xpm*tx+ypm*ty))
                      expdtx=exp(eki*xpm*dtx)
                    else
                      expom=expom*expdtx
                    endif
                  endif

                  wig(ix,iy,itx,ity)=wig(ix,iy,itx,ity)+dreal(em*ep*expom)*wignor

                enddo
              enddo

            enddo
          enddo

        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      end
