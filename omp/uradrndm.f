*CMZ :  3.02/04 13/03/2015  10.24.40  by  Michael Scheer
*CMZ :  3.02/00 28/08/2014  15.45.49  by  Michael Scheer
*CMZ :  2.70/00 12/11/2012  11.53.09  by  Michael Scheer
*CMZ :  2.68/03 29/08/2012  09.57.11  by  Michael Scheer
*-- Author :    Michael Scheer   23/08/2012
      subroutine uradrndm(rn)

*KEEP,gplhint.
*KEND.

c NO WARRANTY

      implicit none
      real rn

      save

      call random_number(rn)

      return
      end
