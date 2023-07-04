*CMZ :  4.00/17 15/11/2022  10.09.59  by  Michael Scheer
*CMZ :  4.00/07 11/05/2020  16.00.17  by  Michael Scheer
*CMZ :  3.07/01 27/03/2019  10.40.38  by  Michael Scheer
*-- Author :    Michael Scheer   27/03/2019
      subroutine waveenvironment

      !use waveenv

      implicit none
*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.

      call get_environment_variable("HOME",chuserhome)
      call get_environment_variable("WAVE",chwavehome)
      call get_environment_variable("PWD",chwavedir)
      call getlog(chuser)

      return
      end
