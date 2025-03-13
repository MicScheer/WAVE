*CMZ :          12/03/2025  10.28.57  by  Michael Scheer
*-- Author :    Michael Scheer   11/03/2025
      integer function isystem(com)

      implicit none

      character(*) com

      call execute_command_line(com,.true.,isystem)

      return
      end
