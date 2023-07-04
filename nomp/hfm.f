*CMZ :  4.00/15 14/03/2022  10.05.01  by  Michael Scheer
*CMZ :  4.00/14 10/02/2022  17.40.28  by  Michael Scheer
*CMZ :  4.00/13 16/12/2021  12.18.27  by  Michael Scheer
*CMZ :  3.02/00 01/09/2014  11.11.30  by  Michael Scheer
*CMZ :  3.01/06 23/06/2014  16.20.53  by  Michael Scheer
*CMZ :  3.01/05 13/06/2014  09.06.46  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  13.12.28  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  09.54.48  by  Michael Scheer
*-- Author :    Michael Scheer   19/01/2012
      subroutine hfm(nid,buffD)

      use clustermod
      use mhbook_mod

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

      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,ntupinfo.
      include 'ntupinfo.cmn'
*KEND.

      double precision buffD(*)
      integer nid,i,k,nvar
      real buff(100),rlow(100),rhigh(100)
      character(80) tit
      character(10) chtag(100)
      character(7) cname,cname2
      character c

      if (icluster.lt.0) then
        if (nid.eq.3601) then
          write(nscr3601,*)buffd(1:36)
        else if (nid.eq.3600) then
          write(nscr3600,*)buffd(1:2)
        else if (nid.eq.3700) then
          write(nscr3700,*)buffd(1:21)
        else if (nid.eq.4600) then
          write(nscr4600,*)buffd(1:5)
        else if (nid.eq.4700) then
          write(nscr4700,*)buffd(1:12)
        else if (nid.eq.30) then
          write(nscr30,*)buffd(1:28)
        else if (nid.eq.7777) then
          write(nscr7777,*)buffd(1:14)
        endif
      endif
      call mh_filln(nid,buffd)
      return

      if (iroottrees.ge.0) then
      endif

      if (iroottrees.eq.0) return

      return
      end
