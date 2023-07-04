*CMZ : 00.00/07 03/08/2010  12.38.35  by  Michael Scheer
*-- Author :    Michael Scheer   02/08/2010
      subroutine util_rfft1b(npoi,xlam0,x,y,nfour,a0,ac,as)

c +PATCH,//UTIL/FOR
c +DECK,util_rfft1b.

c routine is based on fftpack5:
c http://www.scd.ucar.edu/css/software/fftpack5

c input:
c      nfour: number of harmonics
c            product of a view prime numbers.
c      ac,as:  cosine- and sine-like coefficients

c output:
c
c        x,y: function y(x) according to coefficients
c        a0/2+ac(l)*cos(l*k*x)i))+as(i)*sin(l*k*x)
c        x must cover the full intervall lambda=2*pi/k

      implicit none

      double precision x(npoi),y(npoi),ac(npoi),as(npoi),a0,xk,a02,xlam0,
     &  cosn,sinn,cos1,sin1,cosd,sind,xkx
      integer i,nfour,npoi,n

      include 'util_phycon_incl.f'

      xk=2.0d0*pi/xlam0
      a02=a0/2.0d0

      do i=1,npoi
        xkx=xk*x(i)
        cos1=cos(xkx)
        sin1=sin(xkx)
        y(i)=a02+ac(1)*cos1+as(1)*sin1
        cosn=cos1
        sinn=sin1
        do n=2,nfour
          cosd=cosn*cos1-sinn*sin1
          sind=cosn*sin1+sinn*cos1
          sinn=sind
          cosn=cosd
          y(i)=y(i)+ac(n)*cosn+as(n)*sinn
        enddo
      enddo

      return
      end
