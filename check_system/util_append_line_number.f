*CMZ : 00.00/20 22/02/2017  13.35.27  by  Michael Scheer
*-- Author :    Michael Scheer   20/02/2017
      subroutine util_append_line_number

      implicit none

      integer n,luni,luno,last
      character(2048) cline

      last=1
      n=0
      open(newunit=luni,file="util_append_line_number.in")
      open(newunit=luno,file="util_append_line_number.out")
      do while(last.gt.0)
        call util_read_line(luni,cline,last)
        n=n+1
        if (last.gt.0) write(luno,*)cline(1:last)," ",n
      enddo
      close(luno)
      close(luni)

      return
      end
