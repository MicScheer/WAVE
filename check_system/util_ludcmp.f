*CMZ : 00.00/15 24/10/2012  14.26.34  by  Michael Scheer
*CMZ : 00.00/02 27/06/2005  16.08.32  by  Michael Scheer
*-- Author :    Michael Scheer   24/06/2005
      subroutine ludcmp(a,n,np,indx,d,ifail)

c According to Numerical Recipies page 35
c see also util_lubksb

      implicit none

      integer nmax,ifail
      double precision tiny
      parameter (nmax=100,tiny=1.0e-20)

      integer n,np,indx(n),i,j,k,imax
      double precision a(np,np),vv(nmax),d(n),aamax,sum,dum

      stop '*** tut es nicht, (am 24.10.2012 korrigiert!??)***'

      d=1.0
      ifail=0

      if (n.gt.nmax) then
        ifail=1
        goto 999
      endif

      do i=1,n
        aamax=0.
        do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo !j
        if (aamax.eq.0.0d0) then
          ifail=2
          goto 999
        endif
        vv(i)=1./aamax
      enddo !i

      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo !k
          a(i,j)=sum
        enddo !i
        aamax=0.0
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo !k
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        enddo !i
        if (j.ne.imax) then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo !k
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if (a(j,j).eq.0.0) then
          ifail=3
          goto 999
        endif
        if (j.ne.n) then
          dum=1./a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo !i
        endif
      enddo !j

999   return
      end
