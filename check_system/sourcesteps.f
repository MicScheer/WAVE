*CMZ :  3.02/00 10/09/2014  11.06.20  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.66/20 06/07/2011  10.49.51  by  Michael Scheer
*CMZ :  2.66/18 02/12/2010  15.09.50  by  Michael Scheer
*-- Author :    Michael Scheer   02/12/2010
      subroutine sourcesteps(isour,nsteps)
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

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEND.

      implicit none

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEND.

      double precision x1,xendsou,roi(nroip),droix
      integer*8 kzaehl
      integer nsteps,iroi,ir1,ir2,isour,nzaehl

      if (nroi.lt.0) then
        droix=(xendsou-x1)/(nroia-1)
        do iroi=1,nroia
          roix(iroi)=x1+(iroi-1)*droix
          roip(iroi)=1.0d0
        enddo
      endif   !(nroi.lt.0)

      roix(1)=roix(1)-1.0d-6
      roix(nroia)=roix(nroia)+1.0d-6

      do iroi=1,nroia
          roi(iroi)=roix(iroi)
      enddo

      kzaehl=0
      x1=sourceao(1,1,isour)
      xendsou=sourceeo(1,1,isour)
      nzaehl=nlpoio

      ir1=-1
      do iroi=1,nroia
        if (roi(iroi).gt.x1.and.roi(iroi).lt.xendsou.and.ir1.eq.-1) then
          ir1=iroi
          goto 11
        endif
      enddo

11    do iroi=1,nroia
        ir2=iroi
        if (roi(iroi).gt.xendsou) then
          roi(iroi)=xendsou
          ir2=ir2-1
          if (roi(ir2).lt.x1) then
            roi(ir2)=x1
          endif
          goto 12
        endif
      enddo

12    continue

      kzaehl=kzaehl+nzaehl*roip(ir2)*(xendsou-roi(ir2))/(xendsou-x1)

      if (ir1.ne.-1) then

        kzaehl=kzaehl+nzaehl*roip(ir1-1)*(roi(ir1)-x1)/(xendsou-x1)

        do iroi=ir1,ir2-1
          if (roi(iroi).gt.x1.or.roi(iroi)+1.lt.xendsou) then
            kzaehl=kzaehl+nzaehl*roip(iroi)*(roi(iroi+1)-roi(iroi))/(xendsou-x1)
          else if (roi(iroi).gt.x1.or.roi(iroi)+1.lt.xendsou) then
            kzaehl=kzaehl+nzaehl*roip(iroi)*(roi(iroi+1)-roi(iroi))/(xendsou-x1)
          endif
        enddo

      endif

      if (kzaehl.gt.2**30) then
        stop '*** Error in sourcesteps: More then 2**30 steps found for source'
      endif

      nsteps=kzaehl+nbaddp

      return
      end
