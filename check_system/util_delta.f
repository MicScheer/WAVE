*CMZ :          24/08/2018  11.14.04  by  Michael Scheer
*-- Author :    Michael Scheer   24/08/2018
      subroutine util_delta(nline)

      implicit none

      double precision x,y,xo,yo
      integer nline,luni,luno,ieof

      nline=0

      open(newunit=luni,file="util_delta.in",status="old")
      open(newunit=luno,file="util_delta.in",status="new")

      call util_skip_comment_end(luni,ieof)
      if (ieof.ne.0) return
      read(luni,*)x,y
      nline=nline+1
      xo=x
      yo=y

      write(luno,*)x,y,x-xo,y-yo

      do while (.true.)
        call util_skip_comment_end(luni,ieof)
        if (ieof.ne.0) exit
        read(luni,*)x,y
        write(luno,*)x,y,x-xo,y-yo
        nline=nline+1
        xo=x
        yo=y
      enddo

      close(luno)
      close(luni)

      return
      end
