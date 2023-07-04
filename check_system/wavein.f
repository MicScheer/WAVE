*CMZ :  2.16/07 13/09/2000  10.51.26  by  Michael Scheer
*-- Author :    Michael Scheer   21/06/98

c Program to strip WAVE.IN for WAVE under Linux

      program wavein

c reads input file wave.tmp, strips comments etc. and writes wave.in.linux

      implicit none

      integer i,j,ic,icr,iblank

      character*132 c132
      character*1 c1(4)

      equivalence(ic,c1(1))

      data ic/0/

      call system('rm wave.in.linux')
      open(unit=21,file='wave.in.linux',status='new')
      open(unit=20,file='wave.tmp',status='old')

1       read(20,'(a)',end=90)c132

      do i=1,132
        c1(1)=c132(i:i)
        if (ic.ge.65.and.ic.le.90) then
          ic=ic+32
          c132(i:i)=c1(1)
        endif
      enddo

c check, beginning of namelist
c      if (c132(1:1).eq.'$') then
c2       read(20,'(a)',end=90)c132
c        do i=1,132
c             c1(1)=c132(i:i)
c             if (ic.ge.65.and.ic.le.90) then
c               ic=ic+32
c               c132(i:i)=c1(1)
c             endif
c        enddo
c        do i=1,129
c          if (c132(i:i+3).eq.'$end') then
c            goto 1
c          endif
c        enddo
c        goto 2
c      endif   !(c132(1:1).eq.'$')

      icr=0
      do i=1,132
        c1(1)=c132(i:i)
        if (ic.eq.13) then !find carriage return
          icr=i-1
          goto 12
          endif
      enddo

12    backspace(20)

      read(20,'(a)',end=90)c132(1:icr)

      j=icr
      do i=1,icr
        if (c132(i:i).eq.'!') then
          j=i-1
          goto 10
          endif
        enddo

10    continue

      iblank=1
      do i=1,j
        c1(1)=c132(i:i)
        if (c132(i:i).ne.' '.and.ic.ne.9) then
          iblank=0
          goto 20
          endif
        enddo

20    continue

      if (iblank.ne.0) goto 1 !skip blank line

      if (j.lt.icr) then
        write(21,'(a)')c132(1:j)
      else
        write(21,'(a)')c132(1:icr)
      endif

      goto 1   !next line

90    close(20)
      close(21)

      stop
      end
