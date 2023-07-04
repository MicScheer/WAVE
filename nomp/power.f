*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  12.39.05  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  16.46.29  by  Michael Scheer
*CMZ :  2.51/00 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.46/02 07/03/2003  10.27.39  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.19.09  by  Michael Scheer
*CMZ :  2.13/07 17/02/2000  15.11.12  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.01.43  by  Michael Scheer
*CMZ :  2.12/03 07/07/99  12.30.49  by  Michael Scheer
*CMZ :  2.10/01 25/02/99  16.08.38  by  Michael Scheer
*-- Author :    Michael Scheer   25/02/99

      SUBROUTINE POWER(ISOUR)
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
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

        IMPLICIT NONE

        INTEGER ISOUR,ICAL,IOBSV,IOBSVZ

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

        DATA ICAL/0/

        IF (ICAL.EQ.0) THEN

          DO IOBSV=1,NOBSV
            SPECPOWT(IOBSV)=0.D0
          ENDDO

          DO IOBSVZ=1,NOBSVZ
            SPECPOWVT(IOBSVZ)=0.D0
          ENDDO

          SPECPOWVHT=0.D0

          ICAL=1

        ENDIF  !ICAL

        IF (IPIN.NE.0.AND.IPIN.NE.2.and.ipin.ne.3) THEN

          DO IOBSVZ=1,NOBSVZ
            CALL BLENDEPOWV(ISOUR,IOBSVZ)
          ENDDO

          CALL BLENDEPOWVH(ISOUR)

          DO IOBSV=1,NOBSV
            SPECPOWT(IOBSV)=
     &        SPECPOWT(IOBSV)+SPECPOW(ISOUR+NSOURCE*(IOBSV-1))
          ENDDO

          DO IOBSVZ=1,NOBSVZ
            ILIOBZ=ISOUR+NSOURCE*(IOBSVZ-1)
            SPECPOWVT(IOBSVZ)=SPECPOWVT(IOBSVZ)+SPECPOWV(ILIOBZ)
          ENDDO

          SPECPOWVHT=SPECPOWVHT+SPECPOWVH(ISOUR)

        ELSE   !IPIN

          DO IOBSV=1,NOBSV

            SPECPOWT(IOBSV)=SPECPOWT(IOBSV)+
     &        SPECPOW(ISOUR+NSOURCE*(IOBSV-1))

            if (ipin.eq.3) then
              if (ipincirc.eq.0) then
                specpowvh(isour)=specpowvh(isour)+
     &            SPECPOW(ISOUR+NSOURCE*(IOBSV-1))
     &            *pinw*pinh
              else
                specpowvh(isour)=specpowvh(isour)+
     &            SPECPOW(ISOUR+NSOURCE*(IOBSV-1))
     &            *pinr*twopi1
              endif
            endif

          ENDDO   !NOBSV

          SPECPOWVHT=SPECPOWVHT+SPECPOWVH(ISOUR)

        ENDIF  !IPIN

      RETURN
      END
