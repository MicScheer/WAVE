*CMZ :          22/06/2017  12.43.46  by  Michael Scheer
*CMZ : 00.00/07 12/05/2010  13.07.27  by  Michael Scheer
*-- Author :    Michael Scheer   12/04/2010
      subroutine coherentwsum

      use coherentmod

      implicit none

      integer ibunch,iein,ical,icoef

      real, dimension (:), allocatable :: xran

      double complex, dimension (:), allocatable :: field
      double precision, dimension (:), allocatable :: flux

      double complex ampg,x
      double precision phig,phie,wsum,weight

*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      if (ical.eq.0) then

        allocate(xran(nein))

        allocate(flux(nbunch))
        allocate(field(nbunch))

        CALL GRNDMQ(1,1,-1,' ') !S. 39
c        CALL RMARIN(IAMPSEED,NTOTIN,NTOT2IN) !CERN V113

      endif

      flux=0.0d0
      field=(0.0d0,0.0d0)

      do ibunch=1,nbunch

        if (mode.eq.0) then
          call grndm(xran,nein)
        else
          call rnorml(xran,nein)
        endif

        wsum=0.0d0
        do iein=1,nein

          if (mode.eq.0) then
            x=(xran(iein)-0.5)*xlen
            phie=x/xlenfou*twopi1
          else if (mode.eq.0) then
            x=xran(iein)*xlenfou
            phie=x/xlen*twopi1
          else
          endif

          phig=x/xlam*twopi1
          ampg=dcmplx(cos(phig),sin(phig))

          if (mode.eq.0) then
            weight=c(0)
            do icoef=1,ncoef
              weight=weight+c(icoef)*cos(icoef*phie)
            enddo
          else
            weight=1.0d0
          endif

          wsum=wsum+weight
          field(ibunch)=field(ibunch)+weight*ampg

          if (ical.eq.0) then
            if (ibunch.eq.1) then
              write(98,'(5e15.5e3)')
     &          real(ampg),imag(ampg),phie/twopi1,phig/twopi1,weight
            endif
          endif
        enddo !nein

        field(ibunch)=field(ibunch)/sqrt(wsum)
        flux(ibunch)=(dreal(field(ibunch))**2+dimag(field(ibunch))**2)

      enddo !nbunch

      fmean=0.0d0
      frms=0.0d0

      do ibunch=1,nbunch
        fmean=fmean+flux(ibunch)
        frms=frms+flux(ibunch)**2
        if (ical.eq.0) then
          write(99,'(4e15.5e3)')
     &      real(field(ibunch)),imag(field(ibunch)),
     &      flux(ibunch),1./nbunch
        endif
      enddo !nbunch

      fmean=fmean/nbunch
      frms=sqrt(frms/nbunch-fmean**2)

      if (nbunch.eq.1) then
        frms=1.0d0+nein/sqrt(twopi1)
        print*,
     &    '*** Warning: The RMS is estimated from frms=1.0d0+nein/sqrt(2*pi),'
        print*,
     &    'which is only correct, if the spectrum is completely bunched for'
        print*,
     &    'xlamb = xlenfour'
      endif

      ical=1

      return
      end
