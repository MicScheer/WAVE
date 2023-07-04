*CMZ : 00.00/07 04/03/2010  21.48.04  by  Michael Scheer
*-- Author :    Michael Scheer   04/03/2010
      subroutine util_bessel(n,x,bessel,jfail)

c Calculates BESSEL-function Jn(x)

      implicit none

      integer ndimp
      parameter(ndimp=1000)

      double complex z
     &  ,f(ndimp),g(ndimp),fp(ndimp),gp(ndimp),sig(ndimp),eta(ndimp),zlmin

      double precision x,bessel

      integer n,jfail

      if (n.ge.ndimp) then
        jfail=-9999
        return
      endif

      z=dcmplx(x,0.0d0)

c     CALL WCLBES(Z,ETA,ZLMIN,NL,F,G,FP,GP,SIG,KFN,MODE,JFAIL,JPR)

      zlmin=dcmplx(dble(n),0.0d0)
      call wclbes(z,eta,zlmin,0,f,g,fp,gp,sig,2,1,jfail,0) !CERN C309

      bessel=dreal(f(1))

      return
      end
