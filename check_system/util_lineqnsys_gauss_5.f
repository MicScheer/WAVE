*CMZ : 00.00/02 06/11/2006  10.46.35  by  Michael Scheer
*-- Author :    Michael Scheer   24/06/2005
      subroutine util_lineqnsys_gauss_5(a,ndim,b,ws,ifail)

c matrix is diaginal matrix with 5 entries per line, i.e. a(i,i:i+4).ne.0

      integer ndim,ifail,i1,i2,i3,iend

      double precision a(ndim,ndim),b(ndim),ws(ndim),bs

      do i1=ndim,1,-1
        ifail=-1
        if (a(i1,i1).ne.0.0d0) then
          ifail=0
          b(i1)=b(i1)/a(i1,i1)
          iend=min(ndim,i1+4)
          a(i1,i1:iend)=a(i1,i1:iend)/a(i1,i1)
        else
          return
        endif
        do i3=max(1,i1-4),i1-1
c upper triangle
c          iend=1
          b(i3)=b(i3)-b(i1)*a(i3,i1)
c          a(i3,i1:iend)=a(i3,i1:iend)-a(i1,i1:iend)*a(i3,i1)
          a(i3,i1)=a(i3,i1)-a(i1,i1)*a(i3,i1)
        enddo
      enddo

999   return
      end
