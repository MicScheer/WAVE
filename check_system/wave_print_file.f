*CMZ :  4.00/04 17/05/2019  12.18.45  by  Michael Scheer
*-- Author :    Michael Scheer   17/05/2019
      subroutine wave_print_file(luno,chfile)

      implicit none

      integer luno,lunin

      character(*) chfile
      character(1024) cline

      open(newunit=lunin,file=trim(chfile),status='old',err=999)
      do while (.true.)
        read(lunin,'(a)',end=99,err=99) cline
        write(luno,'(a)') trim(cline)
      enddo
99    close(lunin)

999   return
      end
