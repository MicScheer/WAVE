*CMZ :  4.00/15 27/04/2022  08.15.12  by  Michael Scheer
*CMZ :  3.02/03 27/10/2014  12.27.38  by  Michael Scheer
*-- Author :    Michael Scheer   27/10/2014
      subroutine rnorml(rn,n,rr)

      real rn(*),rr(2)
      call util_random_gauss(n,rn,rr)

      return
      end
