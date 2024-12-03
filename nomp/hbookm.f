*CMZ :  4.00/14 22/12/2021  13.41.56  by  Michael Scheer
*CMZ :  4.00/13 16/12/2021  12.18.27  by  Michael Scheer
*CMZ :  3.03/04 27/10/2017  09.47.16  by  Michael Scheer
*CMZ :  3.02/04 12/12/2014  16.03.51  by  Michael Scheer
*CMZ :  3.02/00 01/09/2014  11.11.30  by  Michael Scheer
*CMZ :  3.01/06 23/06/2014  16.20.53  by  Michael Scheer
*CMZ :  3.01/05 12/06/2014  11.35.31  by  Michael Scheer
*CMZ :  3.01/02 25/10/2013  16.22.19  by  Michael Scheer
*CMZ :  2.68/02 04/07/2012  13.35.03  by  Michael Scheer
*CMZ :  2.67/03 09/05/2012  16.32.43  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  12.52.39  by  Michael Scheer
*-- Author :    Michael Scheer   21/01/2012
      subroutine hbookm(nid,ctitle,nvar,cdir,irecl,chtags)
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
      use mhbook_mod
      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,ntupinfo.
      include 'ntupinfo.cmn'
*KEND.

      integer nid,irecl,nvar,i,k,null
      character(7) cname,cname2 ! siehe auch nNTupMax und nHistoMax in waveroot
      character(*) cdir,ctitle
      character(*) chtags(1)
      character c1,cnull

c Wrapper for hbookn, e.g.
c hbookm(NIDTRAC+1,'TRANS. TRACKS',NTUP,'//WAVE',1024,CHTAGS)
c or TTree of root respectively

      equivalence(null,cnull)
      data null/0/

      ! siehe auch nNTupMax und nHistoMax in waveroot
      if (nid.gt.9999)
     & stop '*** Error in HBOOKM: Identifier of Ntuple exceeds 9999'





      call mh_path(cdir)
      call mh_bookn(nid,ctitle,nvar,chtags,irecl)

      return
      end
