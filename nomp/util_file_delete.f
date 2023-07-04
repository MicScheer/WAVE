*CMZ :  3.03/04 02/01/2018  14.56.04  by  Michael Scheer
*CMZ :  3.01/07 23/06/2014  15.51.32  by  Michael Scheer
*CMZ : 00.00/07 21/07/2009  14.58.29  by  Michael Scheer
*CMZ : 00.00/06 12/07/2007  15.45.32  by  Michael Scheer
*-- Author :    Michael Scheer   12/07/2007
      subroutine util_file_delete(file,istat)
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

      integer istat,lun

      character(*) file
      logical lexist,isopen

      istat=-1

      inquire(file=file,exist=lexist)

      if (lexist.eqv..false.) then
        istat=1
        return
      endif

      lun=1234567
      itry=0
1     itry=itry+1
      lun=lun+1
      inquire(unit=lun,opened=isopen)
      if (itry.gt.10) then
        istat=2
        return
      endif
      if (isopen) goto 1

      open(unit=lun,file=file,status='old')
      close(lun,status='delete')

      istat=0

      return
      end
