*CMZ :          26/11/2025  12.10.39  by  Michael Scheer
*-- Author :    Michael Scheer   22/11/2025
*CMZ :  4.01/04 14/11/2023  11.45.54  by  Michael Scheer
*CMZ :  4.01/03 02/06/2023  13.01.26  by  Michael Scheer
*CMZ :  4.00/14 11/02/2022  10.26.56  by  Michael Scheer
*CMZ :  3.05/01 08/05/2018  16.08.32  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  15.45.11  by  Michael Scheer
*CMZ :  2.69/02 07/11/2012  13.59.21  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  12.15.44  by  Michael Scheer
*CMZ :  2.67/00 13/02/2012  10.58.17  by  Michael Scheer
*CMZ :  2.65/03 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.64/04 14/09/2009  15.19.42  by  Michael Scheer
*CMZ :  2.64/03 21/08/2009  17.32.56  by  Michael Scheer
*CMZ :  2.64/02 21/08/2009  17.24.35  by  Michael Scheer
*CMZ :  2.64/01 20/08/2009  15.20.55  by  Michael Scheer
*CMZ :  2.63/05 03/08/2009  16.11.05  by  Michael Scheer
*CMZ :  2.52/12 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.52/11 08/12/2004  13.37.55  by  Michael Scheer
*CMZ :  2.51/02 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.50/00 29/04/2004  15.29.30  by  Michael Scheer
*CMZ :  2.20/01 24/11/2000  21.15.06  by  Michael Scheer
*CMZ :  2.16/08 31/10/2000  14.40.08  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  17.26.23  by  Michael Scheer
*CMZ :  2.13/07 17/02/2000  15.11.12  by  Michael Scheer
*CMZ :  2.13/04 21/01/2000  11.57.55  by  Michael Scheer
*CMZ :  2.13/03 17/01/2000  17.27.08  by  Michael Scheer
*CMZ :  2.12/02 15/06/99  15.16.33  by  Michael Scheer
*CMZ :  2.12/00 27/05/99  10.08.55  by  Michael Scheer
*CMZ :  2.11/01 19/05/99  14.09.50  by  Michael Scheer
*CMZ :  2.10/01 19/03/99  14.13.05  by  Michael Scheer
*CMZ :  2.00/00 06/01/99  11.12.16  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.25.26  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.46.43  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.46  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine reaima_norm

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
*KEEP,reargf90u.
      include 'reargf90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEEP,amplif90u.
      include 'amplif90u.cmn'
*KEND.

      implicit none

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,phasef90.
      include 'phasef90.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEND.

      integer isour,iobsv,kfreq,kmax,ks,kr

      double precision :: smax,rmax,r,specnor_si,s,r1,r2,r3
      double complex :: rea(5),expsh

      if (nsource.gt.1) then
        write(lungfo,*)''
        write(lungfo,*)'        *** Warning in reaima_norm: More then one source points:'
        write(lungfo,*)'        *** Field amplitudes are not correctly calculated due to incoherent ***'
        write(lungfo,*)'        *** treatment of source points, but flux-densities etc. will be ok ***'
        write(lungfo,*)''
        write(6,*)''
        write(6,*)'        *** Warning in reaima_norm: More then one source points:'
        write(6,*)'        *** Field amplitudes are not correctly calculated due to incoherent ***'
        write(6,*)'        *** treatment of source points, but flux-densities etc. will be ok ***'
        write(6,*)''

        return
      endif

      specnor_si= !merke/synchrotron_radiation.txt
     &  dmycur ! strom
     &  /echarge1/hbar1*clight1/PI1*EPS01
     &  *banwid !BW

      expsh=dcmplx(1.0d0,0.0d0)

      do kfreq=1,nfreq

        iobsv=icbrill
        iliobfr=1+nsource*(iobsv-1+nobsv*(kfreq-1))
        iobfr=iobsv+nobsv*(kfreq-1)
        s=spec(iliobfr)
        if (s.gt.0.0d0) then
          if (abs(phgshift).eq.9999.0d0) then
            rea(1:2)=(0.0d0,0.0d0)
            rea(3)=dcmplx(reaima(3,1,iobfr),reaima(3,2,iobfr))
            expsh=rea(3)/abs(rea(3))
            if (phgshift.eq.-9999.0d0) expsh=expsh*cdexp(dcmplx(0.0d0,-pi1/2.0d0))
            do iobsv=1,nobsv
              iobfr=iobsv+nobsv*(kfreq-1)
              rea=dcmplx(reaima(1:5,1,iobfr),reaima(1:5,2,iobfr))/expsh
              reaima(1:5,1,iobfr)=dreal(rea)
              reaima(1:5,2,iobfr)=dimag(rea)
            enddo
          endif
        endif

        do iobsv=1,nobsv
          iliobfr=1+nsource*(iobsv-1+nobsv*(kfreq-1))
          iobfr=iobsv+nobsv*(kfreq-1)
          s=spec(iliobfr)
          if (s.gt.0.0d0) then
            r=
     &        reaima(1,1,iobfr)**2+reaima(1,2,iobfr)**2+
     &        reaima(2,1,iobfr)**2+reaima(2,2,iobfr)**2+
     &        reaima(3,1,iobfr)**2+reaima(3,2,iobfr)**2
            r=sqrt(r/s*specnor_si)
            reaima=reaima/dcmplx(r,1.0d0)
            return
          endif
        enddo
      enddo

      end
