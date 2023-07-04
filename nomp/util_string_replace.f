*CMZ :  4.00/11 04/06/2021  10.29.44  by  Michael Scheer
*-- Author :    Michael Scheer   23/08/2018
      subroutine util_string_replace(cline,chsubstring,chreplace,lensub,lenrep,
     &  n1,n2,istat)

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
c +DECK,util_string_replace.


c Replaces in cline of chsubstring(1:lensub) by chreplace(1:lenrep)
c starting from n1-th to n2-th occurence

      implicit none

      integer istat,l1,l2,nc1,nc2,jstat
      integer lensub,lenrep,leni,leno,n1,n2,ndone

      character(*) cline,chsubstring,chreplace
      !character(2048) cline

      istat=-1

      if (lensub.le.0) then
        leni=len_trim(chsubstring)
      else
        leni=lensub
      endif

      if (lenrep.le.0) then
        leno=len_trim(chreplace)
      else
        leno=lenrep
      endif

      if (len_trim(cline)+leno-leni.gt.len(cline)) return

      jstat=0
      ndone=1

      nc1=max(1,n1)

      if (n2.le.0) then
        nc2=len(cline)
      else
        nc2=n2
      endif

      do ndone=nc1,nc2
        call util_string_substring_counter(cline,chsubstring(1:leni),nc1,
     &    l1,l2,jstat)
        if (jstat.ne.0) exit
        cline=cline(1:l1-1) // chreplace(1:leno)// cline(l2+1:len_trim(cline))
      enddo

      istat=n2-n1+1

      return
      end
