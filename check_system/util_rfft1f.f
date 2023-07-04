*CMZ : 00.00/16 19/03/2014  12.15.48  by  Michael Scheer
*CMZ : 00.00/07 03/08/2010  12.57.29  by  Michael Scheer
*-- Author :    Michael Scheer   02/08/2010
      subroutine util_rfft1f(npoi,y,a0,ac,as,mode,istat)

c +PATCH,//UTIL/FOR
c +DECK,util_rfft1f.

c routine is based on fftpack5:
c http://www.scd.ucar.edu/css/software/fftpack5

c input:
c      npoi: number of data points; routine is more efficient, if npoi is the
c            product of a view prime numbers.
c      y:  data points

C *** X(NPOI)-X(1) IS NOT LAMBDA, BUT LAMBDA = (X(NPOI)-X(1)/(NPOI-1)*NPOI,
C I.E. OPEN INTERVALL!!

c     mode:  -1: new initialization

c output:
c
c        a0: a0/2 is offset
c        ac: cosine-like coefficients
c        as: sine-like coefficients
c
c remark: x,y,a0,ac, and as are double precision variables, but the internal
c         routine is single. If necessary, libfftpack.a might be generated for
c         double precision.
c
c THE OUTPUT-COEFFICIENTS ARE ON CALCULATED FROM 0-npoi/2!!

      implicit none

      double precision y(npoi),ac(npoi),as(npoi),a0

      integer mode,istat,lensav,lenr,lenwrk,inc
      real, dimension (:), allocatable :: wsave,work,r

      integer iallo,i,nh,npoi

      save iallo

      data iallo/0/,inc/1/

      if (mode.eq.-1) then

        if (iallo.ne.0) then
          deallocate (wsave)
          deallocate (work)
          deallocate (r)
          iallo=0
        endif

        if (iallo.eq.0) then
          lensav=npoi + int(log (real(npoi))) + 5
          lenr=inc*(npoi-1)+1
          lenwrk=npoi
          allocate (wsave(lensav))
          allocate (work(lenwrk))
          allocate (r(0:lenr))
          call RFFT1I (npoi, WSAVE, LENSAV, istat)
          if (istat.ne.0) return
          iallo=1
        endif

      endif !mode

      do i=1,npoi
        r(i-1)=y(i)
      enddo

      call rfft1f (npoi, inc, r, lenr, wsave, lensav, work, lenwrk, istat)

      if (istat.ne.0) return

      if (mod(npoi,2).eq.0) then
        nh=npoi/2-1
      else
        nh=(npoi-1)/2
      endif

      a0=2.0d0*r(0)

      do i=1,nh
        ac(i)=r(2*i-1)
        as(i)=r(2*i)
      enddo

      return
      end
