*CMZ : 00.00/19 05/11/2015  15.06.39  by  Michael Scheer
*CMZ : 00.00/15 13/04/2012  12.25.03  by  Michael Scheer
*CMZ : 00.00/11 30/03/2011  14.12.44  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_b_halbach(xperlen,xmaglen,ymaglen,gap,nordmax,br,x,y,b,
     &  eps,
     &  ebeam, harm, defl,
     &  isilent)

      implicit none

      double precision pi,twopi
      parameter (pi=3.14159265359d0)
      parameter (twopi=6.28318530718d0)

      complex*16 z,bc
      double precision xperlen,ymaglen,xmaglen,gap,x,y,b,eps,br,pinm,xk
     &  ,ebeam, harm, defl
      integer nordmax,nu,n,m,isilent

*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,phycon1.
      include 'phycon1.cmn'
*KEND.

      data m/4/

C Calulates magnetic field By for Halbach-Undulator according to
c http://cas.web.cern.ch/cas/Belgium-2009/Lectures/PDFs/Bahrdt-3.pdf
c for M=4

      eps=m*xmaglen/xperlen
      z=dcmplx(x,y)
      xk=twopi/xperlen

      bc=dcmplx(0.0d0,0.0d0)
      do nu=0,nordmax

        n=1+nu*m
        pinm=n*pi/m
        bc=bc
     &    +cos(n*xk*z)
     &    *exp(-n*xk*gap/2.0d0)
     &    *(1.0d0-exp(-n*xk*ymaglen))
     &    *sin(eps*pinm)/pinm

      enddo

      b=2.0d0*br*abs(bc)

      DEFL=ECHARGE1*b*xperlen/1000.0d0/(2.*PI1*EMASSKG1*CLIGHT1)
      HARM=0.95D3*Ebeam**2/(1.D0+DEFL**2/2.D0)/(xperlen/10.0d0)

      if (isilent.eq.0) then
        print*,'Gap, B, K, 1. Harm:',sngl(gap),sngl(b),sngl(defl),sngl(harm)
        print*,'Filling factor:',eps
      endif

      return
      end
