*CMZ :          30/01/2018  13.10.26  by  Michael Scheer
*CMZ : 00.00/19 13/08/2015  13.33.09  by  Michael Scheer
*CMZ : 00.00/07 15/10/2010  09.11.16  by  Michael Scheer
*-- Author :    Michael Scheer   15/10/2010
c +PATCH,//UTIL/FOR
c +DECK,util_getncol.
      subroutine getncol_fortran(file,nskip,ncol)

      implicit none

      integer nskip,ncol
      character*256 file

      open(unit=99,file=file,status='old')
      call getncol(nskip,ncol)
      close(99)

      stop
      end

      subroutine comis_skip_comment_end(lun,ieof)

      implicit none

      integer lun,ieof,icom
      character*10 com,c10

      data c10/'          '/

      icom=0
      ieof=0

1     read(lun,'(a)',end=99) com

      if (
     &  com(1:1).ne.'!'.and.com(1:1).ne.'*'.and.com(1:1).ne.'#'
     & .and.com(1:1).ne.'%'.and.com.ne.''.and.com.ne.c10.and.
     &  com(1:2).ne.' !'.and.com(1:2).ne.' *'.and.com(1:2).ne.' #'
     & .and.com(1:2).ne.' %'
     & ) then
        backspace(lun)
        icom=0
      else
        icom=1
        goto 1
      endif

      return

99    ieof=1

      return
      end


      subroutine getncol(nskip,ncol)

      implicit none

      integer nskip,last,i,ieof,nrow,ncol,inumber,linlen
      CHARACTER*2048 cline
      character c1

      ncol=0

      do nrow=1,nskip
        READ(99,'(A)',END=9)CLINE
      enddo

      CALL comis_skip_comment_end(99,ieof)

      if (ieof.ne.0) return

      READ(99,'(A)',END=9)CLINE

      linlen=len_trim(cline)
      do i=1,linlen
        last=linlen-i
        c1=cline(last:last)
        IF (
     &    C1.EQ.'0'.OR.
     &    C1.EQ.'1'.OR.
     &    C1.EQ.'2'.OR.
     &    C1.EQ.'3'.OR.
     &    C1.EQ.'4'.OR.
     &    C1.EQ.'5'.OR.
     &    C1.EQ.'6'.OR.
     &    C1.EQ.'7'.OR.
     &    C1.EQ.'8'.OR.
     &    C1.EQ.'9'
     &    ) GOTO 19
      enddo

19    CONTINUE

      ncol=0
      inumber=0

      do i=1,last
        c1=cline(i:i)
        IF (
     &      C1.EQ.'0'.OR.
     &      C1.EQ.'1'.OR.
     &      C1.EQ.'2'.OR.
     &      C1.EQ.'3'.OR.
     &      C1.EQ.'4'.OR.
     &      C1.EQ.'5'.OR.
     &      C1.EQ.'6'.OR.
     &      C1.EQ.'7'.OR.
     &      C1.EQ.'8'.OR.
     &      C1.EQ.'9'.OR.
     &      C1.EQ.'.'.OR.
     &      C1.EQ.'+'.OR.
     &      C1.EQ.'-'.OR.
     &      C1.EQ.'D'.OR.
     &      C1.EQ.'E'.OR.
     &      C1.EQ.'d'.OR.
     &      C1.EQ.'e'
     &      ) then
          if (inumber.eq.0) then
            inumber=1
            ncol=ncol+1
          endif
        else
          if (inumber.eq.1) inumber=0
        endif
      enddo

9     return
      end

      subroutine util_getncol(lun,ncol)

      implicit none

      integer lun,last,i,ieof,ncol,inumber,linlen
      CHARACTER*2048 cline
      character c1

      ncol=0

      call util_skip_comment_end(lun,ieof)

      if (ieof.ne.0) return

      READ(lun,'(A)',END=9)CLINE

      linlen=len_trim(cline)
      do i=1,linlen
        last=linlen-i
        c1=cline(last:last)
        IF (
     &    C1.EQ.'0'.OR.
     &    C1.EQ.'1'.OR.
     &    C1.EQ.'2'.OR.
     &    C1.EQ.'3'.OR.
     &    C1.EQ.'4'.OR.
     &    C1.EQ.'5'.OR.
     &    C1.EQ.'6'.OR.
     &    C1.EQ.'7'.OR.
     &    C1.EQ.'8'.OR.
     &    C1.EQ.'9'
     &    ) GOTO 19
      enddo

19    CONTINUE

      ncol=0
      inumber=0

      do i=1,last
        c1=cline(i:i)
        IF (
     &      C1.EQ.'0'.OR.
     &      C1.EQ.'1'.OR.
     &      C1.EQ.'2'.OR.
     &      C1.EQ.'3'.OR.
     &      C1.EQ.'4'.OR.
     &      C1.EQ.'5'.OR.
     &      C1.EQ.'6'.OR.
     &      C1.EQ.'7'.OR.
     &      C1.EQ.'8'.OR.
     &      C1.EQ.'9'.OR.
     &      C1.EQ.'.'.OR.
     &      C1.EQ.'+'.OR.
     &      C1.EQ.'-'.OR.
     &      C1.EQ.'D'.OR.
     &      C1.EQ.'E'.OR.
     &      C1.EQ.'d'.OR.
     &      C1.EQ.'e'
     &      ) then
          if (inumber.eq.0) then
            inumber=1
            ncol=ncol+1
          endif
        else
          if (inumber.eq.1) inumber=0
        endif
      enddo

9     return
      end
