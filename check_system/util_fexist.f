*CMZ :  4.00/04 01/12/2017  12.08.50  by  Michael Scheer
*-- Author :    Michael Scheer   01/12/2017
      function iutil_fexist(filename)

c +PATCH,//UTIL/FOR
c +DECK,util_fexist.

      integer iutil_fexist
      character(*) filename

      logical lexist

      inquire(file=filename,exist=lexist)

      iutil_fexist=0

      if (lexist.eqv..true.) then
        iutil_fexist=1
      endif

      return
      end
