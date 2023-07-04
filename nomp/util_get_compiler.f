*CMZ :  4.00/04 16/05/2019  14.32.08  by  Michael Scheer
*-- Author :    Michael Scheer   16/05/2019
      subroutine util_get_compiler(chcomp,length)

      use ISO_FORTRAN_ENV

      implicit none

      integer length
      character(*) chcomp
      character (len = :), allocatable :: res

      res=compiler_version()
      chcomp=res
      length=len(res)

      end
