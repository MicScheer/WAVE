*CMZ :          16/05/2019  12.37.48  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_guess_platform(chos)

      implicit none

      integer iwin,ianf,iend

      character(*) chos
      character(1024) cshell,csub,copsy,cpath

      call get_environment_variable('OS',copsy)
      call get_environment_variable('SHELL',cshell)
      call get_environment_variable('PATH',cpath)

      chos='UNKNOWN'

      csub="windows"
      call util_lower_case(copsy)
      call util_string_substring(copsy,trim(csub),ianf,iend,iwin)
      print*,iwin

      if (len_trim(cshell).eq.0.and.iwin.eq.0) then
        chos='WINDOWS'
      else if (len_trim(cshell).ne.0.and.iwin.eq.0) then
        chos='CYGWIN'
      else if (len_trim(cshell).ne.0.and.iwin.ne.0) then
        chos='UNIX'
      endif

      return
      end
