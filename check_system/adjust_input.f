*CMZ :  4.00/07 09/01/2020  14.11.19  by  Michael Scheer
*CMZ :  3.03/02 17/12/2015  16.14.05  by  Michael Scheer
*-- Author :    Michael Scheer   11/12/2015
      subroutine adjust_input(nvar_adjust,chvar_adjust,var_adjust)

c Reads wave_adjust.in,
c sets variables chvar_adjust(1:nvar_adjust) to var_adjust(nvar_adjust),
c and writes wave.in

      implicit none

      double precision var_adjust(1000)

      integer lchvar_adjust(1000),loop,nvar_adjust,ivar,max_adjust,iblank,ic1
      integer i,k,ifound,l

      character(128) chvar_adjust(1000)
      character(1024) command_adjust, chfile_adjust,cpwd
      character(512) cline1
      character(256) cline
      character c1
      equivalence(c1,ic1)

      open(unit=99,file='wave.in')
      open(unit=98,file='wave_adjust.in',status='old')

2     read(98,'(a)',end=89)cline

      l=len_trim(cline)

      do i=1,l
        c1=cline(i:i)
        if (c1.eq."!") then
          write(99,'(a)')cline(1:len_trim(cline))
          goto 2
        else if (ic1.eq.32.or.ic1.eq.9) then
          continue
        else
          goto 21
        endif
      enddo

21    do ivar=1,nvar_adjust

        lchvar_adjust(ivar)=len_trim(chvar_adjust(ivar))

        ifound=0

        do i=1,l
          if(cline(i:i+lchvar_adjust(ivar)-1).eq.trim(chvar_adjust(ivar))) then
            ifound=i
            goto 8
          endif
        enddo

8       if (ifound.eq.0.and.ivar.eq.1) then
          write(99,'(a)')cline(1:len_trim(cline))
        else if(ifound.ne.0) then
          do k=i,l
            if (cline(k:k).eq."=") then
              goto 81
            endif
          enddo
81        cline1(1:k)=cline(1:k)
          write(cline1(k+1:k+20),'(g20.14)') var_adjust(ivar)
          cline1(k+21:k+21+l)=cline(i+lchvar_adjust(ivar):l)
          ifound=0
          do i=k,l
            if (cline(i:i).eq."!") then
              ifound=i
              goto 82
            endif
          enddo
82        if (ifound.ne.0) then
            cline1(k+21:k+22)=' '
            cline1(k+22:k+22+l)=cline(i:l)
          endif
          write(99,'(a)')cline1(1:len_trim(cline1))
        endif

      enddo !nvar_adjust
      goto 2

89    close(98)
      close(99)

      return
      end
