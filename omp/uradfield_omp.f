*CMZ :  4.01/04 25/11/2023  13.39.02  by  Michael Scheer
*CMZ :  4.01/02 09/05/2023  13.15.30  by  Michael Scheer
*CMZ :  4.00/15 28/04/2022  11.45.17  by  Michael Scheer
*CMZ :  4.00/13 16/11/2021  17.18.53  by  Michael Scheer
*CMZ :  3.05/05 09/07/2018  15.22.23  by  Michael Scheer
*CMZ :  3.05/04 05/07/2018  08.56.42  by  Michael Scheer
*CMZ :  3.04/00 23/01/2018  17.17.28  by  Michael Scheer
*CMZ :  3.03/04 04/12/2017  15.56.53  by  Michael Scheer
*CMZ :  3.03/02 18/11/2015  12.56.22  by  Michael Scheer
*CMZ :  3.02/04 13/03/2015  10.36.11  by  Michael Scheer
*CMZ :  2.70/00 12/11/2012  11.53.09  by  Michael Scheer
*CMZ :  2.68/04 04/09/2012  09.38.42  by  Michael Scheer
*CMZ :  2.68/03 29/08/2012  12.25.27  by  Michael Scheer
*-- Author :    Michael Scheer   22/08/2012
      subroutine uradfield_omp(x,y,z,bxout,byout,bzout,ex,ey,ez,gamma,istatus,
     &  modewave)

      implicit none

*KEEP,ampli.
      include 'ampli.cmn'
*KEND.

      double precision :: x,y,z,bx,by,bz,ex,ey,ez,
     &  bxout,byout,bzout,gamma,axout,ayout,azout

      integer :: istatus,modewave

      ex=0.0d0
      ey=0.0d0
      ez=0.0d0

      bxout=0.0d0
      byout=0.0d0
      bzout=0.0d0

      if (phrb0v.ne.0.0d0) then

        call bhalba_omp(phrb0v,phrperl,x+phrshift/2.0d0,y,z,bx,by,bz)

        bxout=bxout+bx
        byout=byout+by
        bzout=bzout+bz

      endif

      if (phrb0h.ne.0.0d0) then

        call bhalba_omp(phrb0h,phrperl,x-phrshift/2.0d0,y,z,bx,by,bz)

        bxout=bxout+bx
        byout=byout+bz
        bzout=bzout-by

      endif

      istatus=0

      return
      end
