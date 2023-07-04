*CMZ :  4.00/15 14/03/2022  09.02.26  by  Michael Scheer
*CMZ :  3.02/00 28/08/2014  08.52.10  by  Michael Scheer
*CMZ :  3.01/02 24/02/2014  16.33.26  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.66/17 23/11/2010  09.49.12  by  Michael Scheer
*CMZ :  2.63/05 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.16/08 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.16/07 08/09/2000  17.05.43  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  2.13/03 12/01/2000  11.08.59  by  Michael Scheer
*CMZ :  2.10/01 03/03/99  17.22.37  by  Michael Scheer
*-- Author :    Michael Scheer   24/02/99

      SUBROUTINE PHASE_BAHRDT
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

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,strings.
      include 'strings.cmn'
*KEND.

      INTEGER I,J,K,IFREQ,ifirst,ilast
      DOUBLE PRECISION X,Y,Z,W
      character(64) file,cfreq

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEND.

      DO IFREQ=1,NFREQ

        write(cfreq,*)ifreq
c        ifirst=icfnbl(cfreq,1,len(cfreq))
c        ilast=len_trim(cfreq)
        cfreq=adjustl(cfreq)
        ifirst=1
        ilast=len_trim(cfreq)

        file='eyre-'//cfreq(ifirst:ilast)//'.dat'
        open(unit=99,FILE=file,status='new')

        write(99,*)nobsvz,nobsvy

        do i=1,nobsvy
          do j=1,nobsvz
            k=j+(i-1)*nobsvz
            x=obsv(3,k)
            y=obsv(2,k)
            iobfr=k+nobsv*(ifreq-1)
            z=reaima(2,1,iobfr)
            write(99,*)x*1000.0d0,y*1000.0d0,z
          enddo
        enddo
        close(99)

        file='eyim-'//cfreq(ifirst:ilast)//'.dat'
        open(unit=99,FILE=file,status='new')

        write(99,*)nobsvz,nobsvy
        do i=1,nobsvy
          do j=1,nobsvz
            k=j+(i-1)*nobsvz
            x=obsv(3,k)
            y=obsv(2,k)
            iobfr=k+nobsv*(ifreq-1)
            w=reaima(2,2,iobfr)
            write(99,*)x*1000.0d0,y*1000.0d0,w
          enddo
        enddo
        close(99)

        file='ey2-'//cfreq(ifirst:ilast)//'.dat'
        open(unit=99,FILE=file,status='new')
        write(99,*)nobsvz,nobsvy
        do i=1,nobsvy
          do j=1,nobsvz
            k=j+(i-1)*nobsvz
            iobfr=k+nobsv*(ifreq-1)
            write(99,*)obsv(3,k),obsv(2,k),
     &        reaima(2,1,iobfr)**2+reaima(2,2,iobfr)**2
          enddo
        enddo
        close(99)

        file='ezre-'//cfreq(ifirst:ilast)//'.dat'
        open(unit=99,FILE=file,status='new')
        write(99,*)nobsvz,nobsvy
        do i=1,nobsvy
          do j=1,nobsvz
            k=j+(i-1)*nobsvz
            x=obsv(3,k)
            y=obsv(2,k)
            iobfr=k+nobsv*(ifreq-1)
            z=reaima(3,1,iobfr)
            write(99,*)x*1000.0d0,y*1000.0d0,z
          enddo
        enddo
        close(99)

        file='ezim-'//cfreq(ifirst:ilast)//'.dat'
        open(unit=99,FILE=file,status='new')
        write(99,*)nobsvz,nobsvy
        do i=1,nobsvy
          do j=1,nobsvz
            k=j+(i-1)*nobsvz
            x=obsv(3,k)
            y=obsv(2,k)
            iobfr=k+nobsv*(ifreq-1)
            w=reaima(3,2,iobfr)
            write(99,*)x*1000.0d0,y*1000.0d0,w
          enddo
        enddo
        close(99)

      ENDDO

      RETURN
      END
