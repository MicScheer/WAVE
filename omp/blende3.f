*CMZ :          08/09/2026  09.35.45  by  Michael Scheer
*CMZ :  4.02/01 24/11/2025  18.47.11  by  Michael Scheer
*CMZ :  3.08/01 04/04/2019  12.00.57  by  Michael Scheer
*CMZ :  3.07/00 16/03/2019  15.27.40  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.35/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.34/09 24/09/2001  12.09.00  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  12.20.12  by  Michael Scheer
*CMZ :  2.16/08 24/10/2000  14.08.42  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  16.26.19  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.08.20  by  Michael Scheer
*CMZ :  2.13/04 21/01/2000  12.17.14  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.22.27  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.04.41  by  Michael Scheer
*CMZ : 00.02/04 25/02/97  17.36.07  by  Michael Scheer
*CMZ : 00.02/00 10/12/96  18.09.03  by  Michael Scheer
*CMZ : 00.01/09 01/09/95  13.01.10  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.50.50  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.44  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.06  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE blende3(ISOUR,kfreq)

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

      use circpinmod
      use uradphasemod

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER kfreq,ISOUR,IOBSV,ICAL,lz,ly,nr

      wflux(ISOUR+NSOURCE*(kfreq-1))=0.0d0

      iobsv=0
      nr=0

      do ly=1,nobsvy
        do lz=1,nobsvz

          iobsv=iobsv+1

          if (nrad_u(iobsv).eq.0) cycle

          iobfr=iobsv+nobsv*(kfreq-1)

          if ((
     &        ly.ge.(nobsvy-mobsvy)/2+1.and.ly.le.(nobsvy-mobsvy)/2+mobsvy
     &        .and.
     &        lz.ge.(nobsvz-mobsvz)/2+1.and.lz.le.(nobsvz-mobsvz)/2+mobsvz)
     &        ) then
            nr=nr+nrad_u(iobsv)
            wflux(ISOUR+NSOURCE*(kfreq-1))=wflux(ISOUR+NSOURCE*(kfreq-1))+
     &      spec(isour+nsource*(iobsv-1+nobsv*(kfreq-1)))*pinh*pinw
          endif

          spec(isour+nsource*(iobsv-1+nobsv*(kfreq-1)))=
     &      spec(isour+nsource*(iobsv-1+nobsv*(kfreq-1)))/nrad_u(iobsv)

          spectot(iobfr)=spectot(iobfr)+
     &      spec(isour+nsource*(iobsv-1+nobsv*(kfreq-1)))

        enddo   !iz
      enddo   !iy

      wflux(ISOUR+NSOURCE*(kfreq-1))=wflux(ISOUR+NSOURCE*(kfreq-1))/nr
      wfluxt(kfreq)=wfluxt(kfreq)+wflux(ISOUR+NSOURCE*(kfreq-1))

      RETURN
      END
