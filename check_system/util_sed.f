*CMZ :          05/09/2020  09.15.27  by  Michael Scheer
*-- Author :    Michael Scheer   23/08/2018
      subroutine util_sed(nwords,chwordsin,chwordsout,lenin,lenout,
     &  ncount,chfilein,chfileout,istat)

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

c +PATCH,//UTIL/FOR
c +DECK,util_sed.


c Replaces in lines of chfilein chwordsin(i)(1:lenin(i))
c by chwordsout(i)(1:lenout(i)) starting from ncount(1,i)th to
c ncount(2,i)th occurance and writes result to chfileout
c if ncount(1:2,i)=0 means global replacement

      implicit none

      integer :: isame=0
      integer nwords,lunin,lunout,istat,luntmp,nlines,i,nc1,nc2,jstat
      integer lenin(nwords),lenout(nwords),leni,leno,ncount(2,nwords),ndone

      character(*) chfilein,chfileout
      character(*) chwordsin(nwords),chwordsout(nwords)

      character(2048) cin,cout
      character(128) cwi,cwo

      open(newunit=lunin,file=trim(chfilein),status='old',iostat=istat)
      if (istat.eq.2) return

      istat=0

      if (trim(chfilein).eq.trim(chfileout)) then
        if (nwords.le.0) return
        isame=1
        open(newunit=luntmp,status='scratch')
        nlines=0
        do while (nlines.ge.0)
          read(lunin,'(a)',end=9) cin
          write(luntmp,'(a)') cout
          nlines=nlines+1
        enddo
9       continue
        close(lunin)
        rewind(luntmp)
        lunin=luntmp
      endif !isame

      open(newunit=lunout,file=trim(chfileout),iostat=istat)
      if (istat.ne.2.and.istat.ne.0) goto 9999

      nlines=0
      do while (nlines.ge.0)

        read(lunin,'(a)',end=9999) cin

        cout=cin

        if (nwords.gt.0) then
          do i=1,nwords
            cwi=chwordsin(i)
            cwo=chwordsout(i)
            if (lenin(i).le.0) then
              leni=len_trim(cwi)
            else
              leni=lenin(i)
            endif
            if (lenout(i).le.0) then
              leno=len_trim(cwo)
            else
              leno=lenout(i)
            endif
            jstat=0
            ndone=1
            nc1=max(1,ncount(1,i))
            if (ncount(2,i).le.0) then
              nc2=nwords
            else
              nc2=ncount(2,i)
            endif
            call util_string_replace(cout,cwi(1:leni),cwo(1:leno),leni,leno,
     &        nc1,nc2,jstat)
          enddo !nwords
        endif !(nwords.gt.0) then
        write(lunout,'(a)')trim(cout)
        nlines=nlines+1
      enddo

9999  continue

      if (isame.eq.1) then
        close(luntmp)
      else
        close(lunin)
      endif

      close(lunout)

      return
      end
