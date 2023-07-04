*CMZ : 00.00/02 31/10/2006  09.40.23  by  Michael Scheer
*-- Author :    Michael Scheer   24/06/2005
      subroutine util_lineqnsys_gauss_nmdim(a,ndim,mdim,neqn,b,ws,ifail)

      integer ndim,mdim,neqn,ifail,i1,i2,i3

      double precision a(ndim,mdim),b(mdim),ws(mdim),bs

      do i1=1,neqn
       ifail=-1
       do i2=i1,neqn
          if (a(i2,i1).ne.0.0) then
            bs=b(i2)/a(i2,i1)
          ws(i1:neqn)=a(i2,i1:neqn)/a(i2,i1)
          b(i2)=b(i1)
          a(i2,i1:neqn)=a(i1,i1:neqn)
          b(i1)=bs
          a(i1,i1:neqn)=ws(i1:neqn)
            ifail=0
          goto 9
          endif
         enddo
        goto 999 !return
9     do i3=1,i1-1
              b(i3)=b(i3)-b(i1)*a(i3,i1)
           a(i3,i1:neqn)=a(i3,i1:neqn)-a(i1,i1:neqn)*a(i3,i1)
      enddo
      do i3=i1+1,neqn
           if (a(i3,i1).ne.0.0) then
              b(i3)=b(i3)/a(i3,i1)-b(i1)
           a(i3,i1:neqn)=a(i3,i1:neqn)/a(i3,i1)-a(i1,i1:neqn)
           endif
      enddo
      enddo

999   return
      end
