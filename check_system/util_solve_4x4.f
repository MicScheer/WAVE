*CMZ :          31/10/2022  15.39.11  by  Michael Scheer
*CMZ : 00.00/16 20/07/2015  09.58.29  by  Michael Scheer
*-- Author :    Michael Scheer   20/07/2015
      subroutine util_solve_4x4(a,x,ifail)

      implicit none

      real*8 a(4,4),x(4),det,dws,ws(4,4),xs(4)

      integer ifail,i

      call util_determinante_4(a,det)
      if (det.eq.0.0d0) then
        ifail=-1
        return
      endif

      xs=x

      do i=1,4
        ws=a
        ws(1:4,i)=xs
        call util_determinante_4(ws,dws)
        x(i)=dws/det
      enddo

      ifail=0

      return
      end
