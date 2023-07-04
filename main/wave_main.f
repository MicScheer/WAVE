*CMZ :  4.01/00 13/03/2023  16.34.49  by  Michael Scheer
*CMZ :  4.00/17 15/11/2022  10.50.19  by  Michael Scheer
*CMZ :  4.00/11 30/04/2021  09.31.12  by  Michael Scheer
*CMZ :  4.00/08 07/08/2020  08.27.48  by  Michael Scheer
*CMZ :  4.00/07 18/05/2020  14.17.12  by  Michael Scheer
*CMZ :  4.00/04 05/08/2019  11.45.37  by  Michael Scheer
*CMZ :  3.08/01 02/04/2019  11.52.16  by  Michael Scheer
*CMZ :  3.07/01 21/03/2019  15.16.59  by  Michael Scheer
*CMZ :  3.05/01 04/05/2018  14.52.22  by  Michael Scheer
*CMZ :  3.05/00 26/04/2018  13.14.44  by  Michael Scheer
*CMZ :  3.03/02 29/02/2016  16.23.51  by  Michael Scheer
*CMZ :  3.03/01 05/11/2015  13.09.06  by  Michael Scheer
*CMZ :  3.01/05 13/06/2014  11.46.09  by  Michael Scheer
*CMZ :  3.01/00 17/07/2013  12.16.14  by  Michael Scheer
*CMZ :  3.00/01 20/03/2013  10.21.47  by  Michael Scheer
*CMZ :  3.00/00 14/03/2013  10.40.48  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  15.41.36  by  Michael Scheer
*CMZ :  2.70/00 29/11/2012  16.23.27  by  Michael Scheer
*CMZ :  2.66/07 10/03/2010  16.44.43  by  Michael Scheer
*CMZ :  2.63/02 13/03/2008  15.29.00  by  Michael Scheer
*CMZ :  2.57/05 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.47/23 17/02/2004  13.10.49  by  Michael Scheer
*CMZ :  2.47/09 20/05/2003  14.43.57  by  Michael Scheer
*CMZ :  2.47/08 15/05/2003  16.39.38  by  Michael Scheer
*CMZ :  2.47/03 12/03/2003  15.45.33  by  Michael Scheer
*CMZ :  2.38/00 12/12/2001  16.46.28  by  Michael Scheer
*CMZ :  2.17/00 02/11/2000  16.39.03  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  11.58.00  by  Michael Scheer
*CMZ :  2.13/04 24/01/2000  15.33.01  by  Michael Scheer
*CMZ : 00.01/10 27/08/96  16.19.56  by  Michael Scheer
*CMZ : 00.01/09 09/10/95  17.58.12  by  Michael Scheer
*CMZ : 00.00/07 19/05/94  12.06.35  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  12.01.00  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  17.00.43  by  Michael Scheer

*-- Author : Michael Scheer

*KEEP,GPLHINT.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
*KEND.

c+seq,bpolyederf90m.
c+seq,undumagf90m.
!+seq,waveenv.

      PROGRAM WAVE_MAIN

      use clustermod
      !use waveenv

      implicit none

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEND.
c+seq,platform.
*KEEP,random.
      include 'random.cmn'
*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.

      double precision bound1_adjust(1000),bound2_adjust(1000),
     &  conv_adjust,var_adjust(1000),chimin,chi,
     &  a,b,y,y0,x,x0,xpar(3),ypar(3),apar(3),yp(3),
     &  x1,y1,xopt,yopt,dybad

      integer lchvar_adjust(1000),loop,nvar_adjust,ivar,max_adjust,kvar_adjust,
     &  iblank,ic1,ibad

      integer lenc, nargs
      integer ipid,getpid,i,k,ifound,l,istat

      CHARACTER(8) DAY
      CHARACTER(10) TIM
      character(128) chvar_adjust(1000)
      character(1024) command_adjust, chfile_adjust,cpwd,chrunwave
      character(512) cline1
      character(256) cline

      character c1
      equivalence(c1,ic1)

      CALL ZEIT(6)

      call wavesystem

      chrunwave=trim(chwavehome) // chpathsep  // "bin" //  chpathsep // "wave.exe"

      open(unit=99,file='wave.status')
      write(99,*)"1"
      flush(99)
      close(99)

      ipid=getpid()
      open(unit=99,file='wave.pid')
      write(99,*)ipid
      flush(99)
      close(99)

      nargs=command_argument_count()

      if (nargs.gt.0) then

        iadjust=1

        if (nargs.eq.1) then

          call getcwd(cpwd)
          call get_command_argument(1,chfile_adjust)

          open(unit=99,file=trim(chfile_adjust),status='unknown')
          call util_skip_comment(99)
          read(99,*)nvar_adjust
          if (nvar_adjust.gt.1000) then
            print*,'*** Error in wave_main: Too many parameters to fit in '
            print*,trim(chfile_adjust)
            stop
          endif

          kvar_adjust=nvar_adjust
          if(nvar_adjust.lt.0) kvar_adjust=1

          do ivar=1,kvar_adjust
            call util_skip_comment(99)
            read(99,*)chvar_adjust(ivar),
     &        bound1_adjust(ivar),bound2_adjust(ivar)
            lchvar_adjust(ivar)=len_trim(chvar_adjust(ivar))
          enddo
          call util_skip_comment(99)
          read(99,*)conv_adjust
          call util_skip_comment(99)
          read(99,*)max_adjust
          call util_skip_comment(99)
          read(99,*)command_adjust

          flush(99)
          close(99)

          chimin=1.0d99

          do loop=1,max_adjust

            if (nvar_adjust.eq.-1) then

              if (loop.eq.1) x=bound1_adjust(1)
              if (loop.eq.2) x=bound2_adjust(1)
              if (loop.eq.3) x=(bound1_adjust(1)+bound2_adjust(1))/2.0d0

              var_adjust(1)=x

              call adjust_input(kvar_adjust,chvar_adjust,var_adjust)

              call system(trim(chrunwave))
              call system(trim(command_adjust))

              open(unit=90,file='wave.status',status='old')
              read(90,*) istat
              close(90)

              if (istat.ne.0) then
                print*,"*** Error in adjustment loop: WAVE returned error flag ***"
                stop
              endif

              open(unit=97,file='wave_adjust.dat',status='old')
              if (loop.eq.1) then
                read(97,*) y
                xpar(1)=x
                ypar(1)=y
                open(unit=96,file='wave_adjust.out',status='unknown')
                write(96,*)x,y
                flush(96)
                close(96)
              else if (loop.eq.2) then
                read(97,*) y
                xpar(2)=x
                ypar(2)=y
                open(unit=96,file='wave_adjust.out',status='old',access='append')
                write(96,*)x,y
                flush(96)
                close(96)
              else if (loop.eq.3) then
                read(97,*) y
                xpar(3)=x
                ypar(3)=y
                open(unit=96,file='wave_adjust.out',status='old',access='append')
                write(96,*)x,y
                flush(96)
                close(96)
                istat=0
                call util_parabel(xpar,ypar,apar,yp,xopt,yopt,istat)
                if (istat.ne.0) then
                  print*,'*** Warning in wave_main: util_parabel failed ***'
                  call sleep(3)
                  goto 999
                endif
                x=xopt
                if (xopt.lt.bound1_adjust(1)) then
                  print*,'*** Warning in wave_main: Parabolic optimization failed'
                  print*,'*** Optimum assumed to be lower then first  boundary!'
                  print*,'xopt:',xopt
                  call sleep(3)
                  goto 999
                endif
                if (xopt.gt.bound2_adjust(1)) then
                  print*,'*** Warning in wave_main: Parabolic optimization failed'
                  print*,'*** Optimum assumed to be greater then second  boundary!'
                  print*,'xopt:',xopt
                  call sleep(3)
                  goto 999
                endif
              else
                read(97,*) y
                open(unit=96,file='wave_adjust.out',status='old',access='append')
                write(96,*)x,y
                flush(96)
                close(96)
                ibad=1
                dybad=abs(ypar(1))
                if(abs(ypar(2)).gt.dybad) then
                  ibad=2
                  dybad=abs(ypar(2))
                endif
                if(abs(ypar(3)).gt.dybad) then
                  ibad=3
                endif
                xpar(ibad)=x
                ypar(ibad)=y
                istat=0
                call util_parabel(xpar,ypar,apar,yp,xopt,yopt,istat)
                if (istat.ne.0) then
                  print*,'*** Warning in wave_main: util_parabel failed ***'
                  call sleep(3)
                  goto 999
                endif
                x=xopt
                if (xopt.lt.bound1_adjust(1)) then
                  print*,'*** Warning in wave_main: Parabolic optimization failed'
                  print*,'*** Optimum assumed to be lower then first  boundery!'
                  print*,'xopt:',xopt
                  call sleep(3)
                  goto 999
                endif
                if (xopt.gt.bound2_adjust(1)) then
                  print*,'*** Warning in wave_main: Parabolic optimization failed'
                  print*,'*** Optimum assumed to be greater then second  boundery!'
                  print*,'xopt:',xopt
                  call sleep(3)
                  goto 999
                endif
              endif
              close(97)

              if (abs(y).lt.chimin) then
                chimin=abs(y)
                open(unit=97,file="wave_adjust_opt.in")
                open(unit=96,file="wave.in")
 32             read(96,'(a)',end=42)cline
                write(97,'(a)')trim(cline)
                goto 32
 42             close(96)
                flush(97)
                close(97)
              endif
              if (abs(y).le.conv_adjust) goto 999

            else if (nvar_adjust.gt.0) then

              do ivar=1,kvar_adjust
                var_adjust(ivar)=bound1_adjust(ivar)+
     &            ran(0)*(bound2_adjust(ivar)-bound1_adjust(ivar))
              enddo

              call adjust_input(kvar_adjust,chvar_adjust,var_adjust)
              call system(trim(chrunwave))
              call system(trim(command_adjust));

              open(unit=90,file='wave.status',status='old')
              read(90,*) istat
              close(90)
              if (istat.ne.0) then
                print*,"*** Error in adjustement loop: WAVE returned error flag ***"
                stop
              endif

              open(unit=97,file='wave_adjust.dat',status='old')
              read(97,*) y
              close(97)
              chi=sqrt(y/kvar_adjust)
              if (chi.lt.chimin) then
                chimin=chi
                open(unit=97,file="wave_adjust_opt.in")
                open(unit=96,file="wave.in")
 3              read(96,'(a)',end=4)cline
                write(97,'(a)')trim(cline)
                goto 3
 4              close(96)
                flush(97)
                close(97)
              endif
              if (loop.eq.1) then
                open(unit=96,file='wave_adjust.out',status='unknown')
              else
                open(unit=96,file='wave_adjust.out',status='old',access='append')
              endif
              write(96,*)var_adjust(1:kvar_adjust),chi
              flush(96)
              close(96)

              if (chi.le.conv_adjust.or.loop.eq.max_adjust) then
                goto 999
              endif

            else !nvar_adjust
              print*,'*** Error  in wave_main: Bad mode for parameters to adjust!'
              stop
            endif !nvar_adjust

          enddo

999       CALL ZEIT(6)

          open(unit=96,file="wave_adjust_opt.in")
          open(unit=97,file="wave.in")
 39       read(96,'(a)',end=49)cline
          write(97,'(a)')trim(cline)
          goto 39
 49       close(96)
          flush(97)
          close(97)

          CALL WAVE

          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*)'--- Adjustment loop of  WAVE  completed ---'
          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*)

          stop ' '
        else
          print*,"*** Argument must be the name of a file for the optimization."
          stop '*** Progam WAVE aborted ***'
        endif !nargs.eq.1

      else
        CALL WAVE
      endif !nargs.gt.0

      BMAXGL=DSQRT(DABS(BMAXGL2))

      write(lungfo,*)' '
      write(lungfo,*)'      Number of generated random numbers:', irancalls
      write(lungfo,*)' '

      IF (IABEND.NE.7) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     LAST CHECK:'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'          BMAXGL:',BMAXGL
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'          XBMAXGL:',XBMAXGL
        WRITE(LUNGFO,*)'          YBMAXGL:',YBMAXGL
        WRITE(LUNGFO,*)'          ZBMAXGL:',ZBMAXGL
        WRITE(LUNGFO,*)
      ENDIF ! (IABEND.NE.2)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     CURRENT NUMBER:',ICODE

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'                    --- Program WAVE completed ---'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

      CALL ZEIT(LUNGFO)
      flush(LUNGFO)
      CLOSE (LUNGFO)

      CALL DATE_AND_TIME(DAY,TIM)

      OPEN(UNIT=99,FILE='wave_run.log',status='unknown',access='append')
      WRITE(99,*)-icode,' ',day,' ',tim
      flush(99)
      close(99)

      l=len_trim(code)
      OPEN(UNIT=99,FILE='wave_code.txt',status='unknown')
      WRITE(99,'(a)')'* '//code(1:l)
      flush(99)
      close(99)

      ipid=-1
      open(unit=99,file='wave.pid')
      write(99,*)ipid
      flush(99)
      close(99)

      open(unit=99,file='wave.status')
      write(99,*)"0"
      flush(99)
      close(99)

      if (icluster.ge.0) then
        CALL ZEIT(6)
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)'--- Program WAVE completed ---'
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)
      endif

      END
