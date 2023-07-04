*CMZ :  3.02/03 03/11/2014  10.43.51  by  Michael Scheer
*-- Author :    Michael Scheer   27/10/2014
      function rndm(r)

      implicit none

      real rndm,r,rn(1)

      call util_random(1,rn)

      rndm=rn(1)

      return
      end
