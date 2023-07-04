*CMZ :  2.50/00 23/03/2004  16.37.20  by  Michael Scheer
*CMZ :  2.48/04 16/03/2004  12.51.39  by  Michael Scheer
*CMZ :  2.48/03 03/03/2004  12.49.39  by  Michael Scheer
*CMZ :  2.47/21 03/12/2003  09.13.24  by  Michael Scheer
*CMZ :  2.47/18 27/11/2003  14.41.06  by  Michael Scheer
*-- Author :    Michael Scheer   27/11/2003
      subroutine uncomnamelist
c Program to strip WAVE.IN for WAVE under UNIX

*KEEP,gplhint.
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


c reads input file wave.tmp, strips comments etc. and writes wave.in.linux

      implicit none

      integer i,j,icr,iblank,iline
      byte ic

      character(132) c132
      character c1

      equivalence(ic,c1)

      data ic/0/

      call system('bin/carriagereturn.exe wave.in wave.tmp')
      call system('rm wave.in.linux')

      open(unit=21,file='wave.in.linux',status='new')
      open(unit=20,file='wave.tmp',status='old')

      iline=0
1     read(20,'(a132)',end=90)c132

      iline=iline+1

      do i=1,132
        c1=c132(i:i)
        if (ic.ge.65.and.ic.le.90) then
          ic=ic+32
          c132(i:i)=c1
        endif
      enddo
      do i=1,132
        c1=c132(i:i)
        if (ic.eq.13) then !find carriage return (13)
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
        c1=c132(i:i)
        if (c132(i:i).ne.' '.and.ic.ne.9.) then
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

      return
      end
