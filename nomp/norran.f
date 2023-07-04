*CMZ :  4.00/15 27/04/2022  08.14.32  by  Michael Scheer
*CMZ :  3.02/03 03/11/2014  12.11.22  by  Michael Scheer
*-- Author :    Michael Scheer   27/10/2014
      subroutine norran(rn,rr)

      real rn(1),rr(2)
      call util_random_gauss(1,rn,rr)

      return
      end
