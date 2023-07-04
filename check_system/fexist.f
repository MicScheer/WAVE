*CMZ :  3.03/02 01/03/2016  17.57.48  by  Michael Scheer
*-- Author :    Michael Scheer   01/03/2016
      subroutine fexist(file,iexist)

      implicit none

      integer iexist
      logical lexist
      character(*) file

      iexist=0

      inquire(file=file,exist=lexist)

      if (lexist.eqv..true.) iexist=1

      return
      end
