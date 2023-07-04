*CMZ :          26/02/2021  13.36.16  by  Michael Scheer
*CMZ : 00.00/03 23/11/2006  10.24.13  by  Michael Scheer
*CMZ : 00.00/02 27/06/2005  14.10.03  by  Michael Scheer
*-- Author :    Michael Scheer   24/06/2005
      subroutine util_lineqnsys_gauss(a,n,b,ws,ifail)

c WARNING: Matrix a will be destroyed

      integer n,ifail,i1,i2,i3

      double precision a(n,n),b(n),ws(n),bs,eps

      data eps/1.0d-20/

      do i1=1,n

        ifail=-1

        do i2=i1,n
         if (abs(a(i2,i1)).ge.eps) then
           bs=b(i2)/a(i2,i1)
           ws(i1:n)=a(i2,i1:n)/a(i2,i1)
           b(i2)=b(i1)
           a(i2,i1:n)=a(i1,i1:n)
           b(i1)=bs
           a(i1,i1:n)=ws(i1:n)
           ifail=0
           goto 9
         endif
       enddo

       return

9      do i3=1,i1-1
         b(i3)=b(i3)-b(i1)*a(i3,i1)
         a(i3,i1:n)=a(i3,i1:n)-a(i1,i1:n)*a(i3,i1)
       enddo
       do i3=i1+1,n
         if (abs(a(i3,i1)).ge.eps) then
           b(i3)=b(i3)/a(i3,i1)-b(i1)
           a(i3,i1:n)=a(i3,i1:n)/a(i3,i1)-a(i1,i1:n)
         endif
       enddo
      enddo

      return
      end
