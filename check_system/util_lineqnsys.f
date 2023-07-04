*CMZ : 00.00/02 27/06/2005  16.08.32  by  Michael Scheer
*-- Author :    Michael Scheer   24/06/2005
      subroutine util_lineqnsys(a,n,ndim,b,ifail)

      integer ndimp
      parameter (ndimp=100)
      integer indx(ndimp),n,ndim,ifail

      double precision a(ndim,ndim),b(n)

      stop '*** tut es nicht ***'

      ifail=0

      if (n.gt.ndimp.or.ndim.gt.ndimp) then
        ifail=-1
        goto 999
      endif

      call ludcmp(a,n,ndim,indx,b,ifail)
      if (ifail.ne.0) goto 999
      call lubksb(a,n,ndim,indx,b)

999   return
      end
