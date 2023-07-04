*CMZ : 00.00/07 10/06/2008  13.28.49  by  Michael Scheer
*CMZ : 00.00/02 07/01/98  15.23.41  by  Michael Scheer
*-- Author :    Michael Scheer   07/01/98
      subroutine util_hbook_ini

      implicit none

      integer ndpawcp
      parameter (ndpawcp=200000)
      real*4 rpaw(ndpawcp)
      common/pawc/rpaw

      call hlimit(ndpawcp)

      return
      end
