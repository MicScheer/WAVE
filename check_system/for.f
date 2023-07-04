*CMZ : 00.00/15 24/10/2012  14.33.54  by  Michael Scheer
*CMZ : 00.00/06 31/10/2007  17.15.57  by  Michael Scheer
*-- Author :    Michael Scheer   30/10/2007
      subroutine util_rfft(n,f,cfft,ifail)

c +PATCH,//UTIL/FOR
c +DECK,util_rfft.

c fft of f(0:n-1), n is power of 2
c
c output: fft(odd)  = coef_cos
c         fft(even) = coef_sin
c         f = fft(1) + fft(2)*cos(k) +fft(4)*cos(2*k) ...
c           +          fft(3)*sin(k) +fft(5)*sin(2*k) ...

      implicit none

      real f(0:n-1)
      complex cfft(0:n/2-1)

      integer n,ifail,ipow,i

      ifail=0
      ipow=nint(alog(float(n))/alog(2.0))

      if (2**ipow.gt.n) ipow=ipow-1

      if (n.ne.2**ipow) then
        ifail=-1
        return
      endif

      do i=0,n/2-1
        cfft(i)=cmplx(f(2*i),f(2*i+1))
      enddo

      call rfstft(-ipow,cfft)

      do i=1,n/2-1
        cfft(i)=cmplx(
     &     real(cfft(i)+conjg(cfft(i))),
     &    -imag(cfft(i)-conjg(cfft(i)))
     &    )
      enddo

      return
      end

