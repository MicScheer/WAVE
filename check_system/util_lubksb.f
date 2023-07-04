*CMZ : 00.00/02 27/06/2005  16.08.32  by  Michael Scheer
*-- Author :    Michael Scheer   24/06/2005
      subroutine lubksb(a,n,np,indx,b)

c According to numerical recipies page 36
c see also util_ludcmp

      implicit none

      integer n,np,indx(n),ii,ll,j,i
      double precision a(np,np),b(n),sum

      stop '*** tut es nicht ***'

      ii=0

      do i=1,n
        ll=indx(n)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0) then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo !j
        else if (sum.ne.0.0) then !ii
          ii=i
        endif !ii
        b(i)=sum
      enddo !i

      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo !j
        b(i)=sum/a(i,i)
      enddo

      return
      end
