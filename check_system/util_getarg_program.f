*CMZ :          04/09/2020  08.27.25  by  Michael Scheer
*-- Author :    Michael Scheer   04/09/2020
      subroutine util_getarg_program(nlen,cprog)

      !Liefert den Namen des laufenden Programmes

      implicit none

      integer nlen,ierr
      character(512) cprog

      cprog=''
      call get_command_argument(0,cprog,nlen,ierr)

      return
      end
