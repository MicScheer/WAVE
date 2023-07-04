*CMZ :  4.00/11 16/05/2021  13.16.14  by  Michael Scheer
*CMZ :  3.03/02 15/12/2015  15.57.23  by  Michael Scheer
*CMZ :  3.01/02 24/01/2014  17.44.47  by  Michael Scheer
*CMZ :  3.01/00 04/07/2013  08.33.03  by  Michael Scheer
*CMZ :  3.00/01 02/04/2013  13.47.23  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.63/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.61/06 12/04/2007  09.47.01  by  Michael Scheer
*CMZ :  2.58/01 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.57/05 24/08/2006  16.44.25  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.35.17  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  14.47.19  by  Michael Scheer
*CMZ : 00.00/07 01/06/94  10.41.53  by  Michael Scheer
*-- Author : Michael Scheer
C-----------------------------------------------------------
      subroutine bfoursincos(xin,yin,zin,
     &  bxout,byout,bzout,axout,ayout,azout,
     &  nfour,xl0,zl0,cfft,init)

C--------------------------------------------------------------------
      ! NOTE: The interpretation and as in util_rfft, NOT AS IN BFOUR
      ! xin is longitudinal
C--------------------------------------------------------------------

      implicit none

      integer :: nfour,iallo=0,nfourold=-1,init,kinit,i,k

      double precision :: xin,yin,zin,bxout,byout,bzout,axout,ayout,azout,
     &  xl0,zl0,bxh,byh,bzh,axh,ayh,azh,a0,dchyky,dcsxkx,dcszkz,dexpomy,
     &  dshyky,dsnxkx,dsnzkz,expomy,expomy1,xk0four,zk0four,
     &  dnull=0.0d0

      double precision, dimension(:), allocatable :: ac,as,xkfour,ykfour,zkfour

      double complex cdexpomx,cdexpomz,cexpomz

      complex cfft(0:nfour-1)

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      save

      if (init.ne.0.or.nfourold.ne.nfour) then

        nfourold=nfour

        if (iallo.ne.0) then
          deallocate(ac,as,xkfour,ykfour,zkfour)
        endif

        allocate(ac(nfour-1),as(nfour-1),
     &    xkfour(nfour-1),ykfour(nfour-1),zkfour(nfour-1))
        iallo=1

        if (xl0.ne.0.0d0) then
          zk0four=twopi1/xl0
        else
          zk0four=0.0d0
        endif

        if (zl0.ne.0.0d0) then
          xk0four=twopi1/zl0
        else
          xk0four=0.0d0
        endif

        do i=1,nfour-1
          zkfour(i)=zk0four*dble(i)
          xkfour(i)=xk0four
          ykfour(i)=dsqrt(zkfour(i)**2+xkfour(i)**2)
        end do

      endif

      a0=dble(real(cfft(0)))
      ac(1:nfour-1)=dble(real(cfft(1:nfour-1)))
      as(1:nfour-1)=dble(imag(cfft(1:nfour-1)))

c if changed, consider following loop and sr mybfeld



      axh= a0/2.0d0*xin !xin is here z
      ayh=0.
      azh=0.

      cdexpomx=cdexp(dcmplx(dnull,xkfour(1)*(-zin)))
      dcsxkx=dreal(cdexpomx)
      dsnxkx=dimag(cdexpomx)

      dexpomy=dexp(ykfour(1)*yin)
      expomy=1.0d0

      cdexpomz=cdexp(dcmplx(dnull,zkfour(1)*xin))
      cexpomz=dcmplx(1.0d0,dnull)

      bxh=0.0d0
      byh=a0
      bzh=0.0d0

      do k=1,nfour-1

        if (xk0four.ne.0.0d0) then
          expomy=dexp(ykfour(k)*yin)
        else
          expomy=expomy*dexpomy
        endif

        expomy1=1.0d0/expomy
        dchyky=(expomy+expomy1)*0.5d0
        dshyky=(expomy-expomy1)*0.5d0

        cexpomz=cexpomz*cdexpomz
        dcszkz=dreal(cexpomz)
        dsnzkz=dimag(cexpomz)

        bxh=bxh-xkfour(k)/ykfour(k)*dsnxkx*dshyky*(ac(k)*dcszkz+as(k)*dsnzkz)
        byh=byh+                    dcsxkx*dchyky*(ac(k)*dcszkz+as(k)*dsnzkz)
        bzh=bzh-zkfour(k)/ykfour(k)*dcsxkx*dshyky*(ac(k)*dsnzkz-as(k)*dcszkz)



        axh=axh+ac(k)/zkfour(k)*dcsxkx*dchyky*dsnzkz
        azh=azh+0.0

        ayh=ayh+ac(k)/zkfour(k)*xkfour(k)/ykfour(k)*dsnxkx*dshyky*dsnzkz

      enddo

      bzout=-bxh
      byout= byh
      bxout= bzh

      azout=-axh
      ayout= ayh
      axout= azh

      return
      end
