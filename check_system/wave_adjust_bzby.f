*CMZ :  3.03/02 11/01/2016  09.50.20  by  Michael Scheer
*-- Author :    Michael Scheer   11/01/2016
      program wave_adjust

c     Utility to write the current results of wave.adjust to wave_adjust.dat

      implicit none

      double precision chi2,shift,by,bz
      character(1024) cline

      open(unit=20,file='wave.out')
1     read(20,'(a)') cline

        if (cline(1:15).eq."       URSHIFT:")
     &    read(cline(16:len_trim(cline)),*) shift

        if (cline(1:18).eq."      BYmax,BYmin:")
     &    read(cline(19:len_trim(cline)),*) by

        if (cline(1:18).eq."      BZmax,BZmin:") then
          read(cline(19:len_trim(cline)),*) bz
          goto 9
        endif

      goto 1

9     close(20)

      chi2=(1.0d0-abs(bz/by))**2

      open(unit=20,file='wave_adjust.dat')
      write(20,*)chi2
      close(20)

      end
