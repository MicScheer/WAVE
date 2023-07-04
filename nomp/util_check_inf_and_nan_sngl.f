*CMZ :  3.05/12 23/08/2018  11.56.03  by  Michael Scheer
*-- Author :    Michael Scheer   23/08/2018
      subroutine util_check_inf_and_nan_sngl(n,x,infnan,istatus)

      implicit none

      real x(*)
      integer n,i,istatus,infnan(*)

      istatus=0

      do i=1,n
        if (x(i).ne.x(i)) then
          infnan(i)=1
          istatus=istatus+1
        else if (1.0/x(i).eq.0.0) then
          infnan(i)=-1
          istatus=istatus+1
        endif
      enddo

      return
      end
