*CMZ :          08/09/2026  13.36.33  by  Michael Scheer
*CMZ :  4.02/01 11/12/2025  22.03.18  by  Michael Scheer
*CMZ :  4.01/05 18/04/2024  13.59.40  by  Michael Scheer
*CMZ :  4.00/15 07/04/2022  07.14.03  by  Michael Scheer
*CMZ :  4.00/14 30/12/2021  15.41.22  by  Michael Scheer
*CMZ :  4.00/13 07/12/2021  18.47.10  by  Michael Scheer
*CMZ :  3.03/02 07/12/2015  17.18.10  by  Michael Scheer
*CMZ :  3.02/06 17/04/2015  16.27.01  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.13.36  by  Michael Scheer
*CMZ :  2.70/11 18/02/2013  16.49.40  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  12.17.08  by  Michael Scheer
*CMZ :  2.68/01 29/05/2012  16.50.03  by  Michael Scheer
*-- Author :    Michael Scheer   29/05/2012
      subroutine hfreq3
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

      use uradphasemod

      implicit none

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,specdip.
      include 'specdip.cmn'
*KEEP,phasef90.
      include 'phasef90.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      COMPLEX*16 e(3),b(3)

      DOUBLE PRECISION WEIGHT,smax,specnor_si,reanor,rn(3),
     &  dist,dist0,ddist,h2,censoux,censouy,censouz,dphase,wlen,waves

      real*8 fstuple(8),fspec(34)

      real df,flow,fhig,ff
      integer ifreq,id,icycle,mfreq,iobsv,iz,iy,isour,iobsvz,iobsvy,iobsvr,iobsvphi,ifrq

      CHARACTER(4) CHSTOK(12)

      if (ipin.ne.3) return

      SPECNOR_SI= !merke/synchrotron_radiation.txt
     &  dmycur ! Strom
     &  /echarge1/hbar1*clight1/PI1*EPS01
     &  *banwid !BW

c      reanor=sqrt(smax/reanor/specnor_si)
      reanor=1.0d0

      DO ISOUR=1,NSOURCE

        if (ispecdip.le.0) then
          censoux=(min(sourceeo(1,1,isour),xiend)
     &      +max(sourceao(1,1,isour),xianf))/2.d0
          censouy=sourceeo(2,1,isour)
          censouz=sourceeo(3,1,isour)
        else
          censoux=x0dip(isour)
          censouy=y0dip(isour)
          censouz=z0dip(isour)
        endif

        dist0=sqrt(
     &    ((pincen(3)-censouz)**2+
     &    (pincen(2)-censouy)**2)+
     &    (pincen(1)-censoux)**2)

        DO IOBSV=1,NOBSV

          IF (IPIN.NE.0) THEN
            IOBSVY=(IOBSV-1)/NOBSVZ+1
            IOBSVZ=IOBSV-NOBSVZ*(IOBSVY-1)
          ELSE
            IOBSVY=0
            IOBSVZ=0
          ENDIF

          DO ifrq=1,NFREQ,IHFREQ
            FSPEC(1)=ISOUR
            FSPEC(2)=IOBSV
            FSPEC(3)=OBSV(1,IOBSV)
            FSPEC(4)=OBSV(2,IOBSV)
            FSPEC(5)=OBSV(3,IOBSV)
            if (abs(fspec(4)).lt.1.0d-15) fspec(4)=0.0d0
            if (abs(fspec(5)).lt.1.0d-15) fspec(5)=0.0d0
            FSPEC(6)=FREQ(ifrq)
            FSPEC(7)=SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(ifrq-1)))
            FSPEC(8)=IOBSVZ
            FSPEC(9)=IOBSVY
            FSPEC(10)=ifrq
            IOBFR=IOBSV+NOBSV*(ifrq-1)
c            IF (ISPECMODE.EQ.3) THEN

            FSPEC(11)=reanor*reaIMA(1,1,IOBFR)
            FSPEC(12)=reanor*reaIMA(1,2,IOBFR)
            FSPEC(13)=reanor*reaIMA(2,1,IOBFR)
            FSPEC(14)=reanor*reaIMA(2,2,IOBFR)
            FSPEC(15)=reanor*reaIMA(3,1,IOBFR)
            FSPEC(16)=reanor*reaIMA(3,2,IOBFR)

            FSPEC(17)=reanor*reaIMA(4,1,IOBFR)
            FSPEC(18)=reanor*reaIMA(4,2,IOBFR)
            FSPEC(19)=reanor*reaIMA(5,1,IOBFR)
            FSPEC(20)=reanor*reaIMA(5,2,IOBFR)

            e(1)=dcmplx(fspec(11),fspec(12))
            e(2)=dcmplx(fspec(13),fspec(14))
            e(3)=dcmplx(fspec(15),fspec(16))

            FSPEC(21)=reanor*reaIMA(6,1,IOBFR)
            FSPEC(22)=reanor*reaIMA(6,2,IOBFR)
            FSPEC(23)=reanor*reaIMA(7,1,IOBFR)
            FSPEC(24)=reanor*reaIMA(7,2,IOBFR)
            FSPEC(25)=reanor*reaIMA(8,1,IOBFR)
            FSPEC(26)=reanor*reaIMA(8,2,IOBFR)

            b(1)=dcmplx(fspec(21),fspec(22))
            b(2)=dcmplx(fspec(23),fspec(24))
            b(3)=dcmplx(fspec(25),fspec(26))

            FSPEC(27)=reanor*reaIMA(9,1,IOBFR)
            FSPEC(28)=reanor*reaIMA(9,2,IOBFR)
            FSPEC(29)=reanor*reaIMA(10,1,IOBFR)
            FSPEC(30)=reanor*reaIMA(10,2,IOBFR)

            rn(1)=real(e(2)*conjg(b(3))-e(3)*conjg(b(2)))
            rn(2)=real(e(3)*conjg(b(1))-e(1)*conjg(b(3)))
            rn(3)=real(e(1)*conjg(b(2))-e(2)*conjg(b(1)))

            rn=rn/norm2(rn)
            fspec(32:34)=rn(1:3)

            dist=sqrt(
     &          ((obsv(3,iobsv)-censouz)**2+
     &        (obsv(2,iobsv)-censouy)**2)+
     &        (obsv(1,iobsv)-censoux)**2)

            ddist=dist-dist0

            wlen=clight1*hbarev1*twopi1/freq(ifrq)
            waves=ddist/wlen
            dphase=waves*twopi1

            FSPEC(31)=dphase

            CALL hfm(NIDSPEC,FSPEC)

          ENDDO   !NFREQ
        ENDDO   !IOBSV
      ENDDO   !ISOUR

      IF (MPINR.NE.0) THEN
        DO ISOUR=1,NSOURCE
          DO IOBSV=1,NOBSVRPHI
            IF (IPIN.NE.0) THEN
              IOBSVPHI=(IOBSV-1)/NOBSVR+1
              IOBSVR=IOBSV-NOBSVR*(IOBSVPHI-1)
            ELSE
              IOBSVR=0
              IOBSVPHI=0
            ENDIF
            DO ifrq=1,NFREQ,IHFREQ
              FSPEC(1)=ISOUR
              FSPEC(2)=IOBSV
              FSPEC(3)=OBSVRPHI(1,IOBSV)
              FSPEC(4)=OBSVRPHI(2,IOBSV)*SIN(OBSVRPHI(3,IOBSV))
              FSPEC(5)=OBSVRPHI(2,IOBSV)*COS(OBSVRPHI(3,IOBSV))
              FSPEC(6)=OBSVRPHI(2,IOBSV)
              FSPEC(7)=OBSVRPHI(3,IOBSV)
              FSPEC(8)=FREQ(ifrq)
              FSPEC(9)=SPECRPHI(ISOUR+NSOURCE*(IOBSV-1+NOBSVRPHI*(ifrq-1)))
              FSPEC(10)=IOBSVR
              FSPEC(11)=IOBSVPHI
              FSPEC(12)=ifrq
              IOBFR=IOBSV+NOBSVRPHI*(ifrq-1)
              FSPEC(13)=reanor*reaIMARPHI(1,1,IOBFR)
              FSPEC(14)=reanor*reaIMARPHI(1,2,IOBFR)
              FSPEC(15)=reanor*reaIMARPHI(2,1,IOBFR)
              FSPEC(16)=reanor*reaIMARPHI(2,2,IOBFR)
              FSPEC(17)=reanor*reaIMARPHI(3,1,IOBFR)
              FSPEC(18)=reanor*reaIMARPHI(3,2,IOBFR)
              FSPEC(19)=reanor*reaIMARPHI(4,1,IOBFR)
              FSPEC(20)=reanor*reaIMARPHI(4,2,IOBFR)
              FSPEC(21)=reanor*reaIMARPHI(5,1,IOBFR)
              FSPEC(22)=reanor*reaIMARPHI(5,2,IOBFR)

              FSPEC(23)=reanor*reaIMARPHI(6,1,IOBFR)
              FSPEC(24)=reanor*reaIMARPHI(6,2,IOBFR)
              FSPEC(25)=reanor*reaIMARPHI(7,1,IOBFR)
              FSPEC(26)=reanor*reaIMARPHI(7,2,IOBFR)
              FSPEC(27)=reanor*reaIMARPHI(8,1,IOBFR)
              FSPEC(28)=reanor*reaIMARPHI(8,2,IOBFR)
              FSPEC(29)=reanor*reaIMARPHI(9,1,IOBFR)
              FSPEC(30)=reanor*reaIMARPHI(9,2,IOBFR)
              FSPEC(31)=reanor*reaIMARPHI(10,1,IOBFR)
              FSPEC(32)=reanor*reaIMARPHI(10,2,IOBFR)

              h2=(obsvr(iobsvr)/dist0)**2

              if (h2.lt.0.01) then
                ddist=dist0*(h2/2.0d0-h2**2/8.0d0)
              else
                ddist=dist0*(sqrt(1.0d0+h2)-1.0d0)
              endif

              dphase=ddist/freq(ifrq)*wtoe1*1.0d9*twopi1

              FSPEC(33)=dphase

              CALL hfm(NIDSPECRPHI,FSPEC)

            ENDDO   !NFREQ
          ENDDO   !IOBSV
        ENDDO   !ISOUR
      ENDIF !MPINR

      do isour=1,nsource
        iobsv=0
        do iy=1,nobsvy
          do iz=1,nobsvz
            iobsv=iobsv+1
            fstuple(1:3)=obsv(1:3,iobsv)
            fstuple(4)=specpow(iobsv+nobsv*(isour-1))
            fstuple(5)=iz
            fstuple(6)=iy
            fstuple(7)=isour
            call hfm(nidpow,fstuple)
          enddo
        enddo
      enddo

      df=freq(2)-freq(1)
      flow=freq(1)-df/2.
      fhig=freq(nfreq)+df/2.

      if (ifreq2p.eq.1.or.freqlow.eq.freqhig) then
          DF=freqhig-freqlow
          FLOW=freqlow-DF/2.
          FHIG=freqlow+DF/2.
      else if (ifreq2p.eq.-1) then
        DF=freqhig-freqlow
        ff=(freqlow+freqhig)/2.
        FLOW=ff-DF/2.
        FHIG=ff+DF/2.
      endif

      if (flow.lt.0.) then
        write(lungfo,*)
        write(lungfo,*)'*** WARNING IN HFREQ3 ***'
        write(lungfo,*)'LOW EDGE OF HISTOGRAM NEGATIVE'
        write(lungfo,*)'BE CAREFUL IF X-AXIS IS PLOTTED WITH LOGARITHMIC SCALE'
        write(lungfo,*)
      endif

      id=icfreq
      mfreq=nint((fhig-flow)/df)
      call hbook1m(id,'Mean flux-density in pinhole x 1.e-6',
     &  mfreq,flow,fhig,vmx)
      do ifrq=1,nfreq,ihfreq
        do iobsv=1,nobsv
          iobfr=iobsv+nobsv*(ifrq-1)
          call hfillm(id,sngl(freq(ifrq)),0.,spectot(iobfr)*1.0d-6)
        enddo
      enddo   !nfreq
      call mhrout(id,icycle,' ')

      id=idfreq
      mfreq=nint((fhig-flow)/df)
      call hbook1m(id,'Photon flux through pinhole',
     &  mfreq,flow,fhig,vmx)
      do ifrq=1,nfreq,ihfreq
        call hfillm(id,sngl(freq(ifrq)),0.,wfluxt(ifrq))
      enddo   !nfreq
      call mhrout(id,icycle,' ')

      do ifrq=1,nfreq
        fstuple(1)=freq(ifrq)
        fstuple(2)=wfluxt(ifrq)
        call hfm(nidfreqp,fstuple)
      enddo   !nfreq
      call mhrout(nidfreqp,icycle,' ')

      if (istokes.ne.0) then

        id=icfrs0
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Mean flux-density S0 through pinhole x 1.e-6',
     &    mfreq,flow,fhig,vmx)
        do ifrq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifrq)),0.,dble(stokes(1,ifrq))*1.0d-6)
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        id=icfrs1
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Mean flux-density S1 through pinhole x 1.e-6',
     &    mfreq,flow,fhig,vmx)
        do ifrq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifrq)),0.,dble(stokes(2,ifrq))*1.0d-6)
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        call mhrout(id,icycle,' ')
        id=icfrs2
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Mean flux-density S2 through pinhole x 1.e-6',
     &    mfreq,flow,fhig,vmx)
        do ifrq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifrq)),0.,dble(stokes(3,ifrq))*1.0d-6)
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        call mhrout(id,icycle,' ')
        id=icfrs3
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Mean flux-density S3 through pinhole x 1.e-6',
     &    mfreq,flow,fhig,vmx)
        do ifrq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifrq)),0.,dble(stokes(4,ifrq))*1.0d-6)
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        id=idfrs0
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Flux S0 through pinhole',
     &    mfreq,flow,fhig,vmx)
        do ifrq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifrq)),0.,dble(wstokes(1,ifrq)))
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        id=idfrs1
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Flux S1 through pinhole',
     &    mfreq,flow,fhig,vmx)
        do ifrq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifrq)),0.,dble(wstokes(2,ifrq)))
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        call mhrout(id,icycle,' ')
        id=idfrs2
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Flux S2 through pinhole',
     &    mfreq,flow,fhig,vmx)
        do ifrq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifrq)),0.,dble(wstokes(3,ifrq)))
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        call mhrout(id,icycle,' ')
        id=idfrs3
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Flux S3 through pinhole',
     &    mfreq,flow,fhig,vmx)
        do ifrq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifrq)),0.,dble(wstokes(4,ifrq)))
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        id=4600
        do ifrq=1,nfreq
          fstuple(1)=freq(ifrq)
          fstuple(2)=wstokes(1,ifrq)
          fstuple(3)=wstokes(2,ifrq)
          fstuple(4)=wstokes(3,ifrq)
          fstuple(5)=wstokes(4,ifrq)
          call hfm(id,fstuple)
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        DO IOBSV=1,NOBSV
          IF (IPIN.NE.0) THEN
            IOBSVY=(IOBSV-1)/NOBSVZ+1
            IOBSVZ=IOBSV-NOBSVZ*(IOBSVY-1)
          ELSE
            IOBSVY=1
            IOBSVZ=1
          ENDIF
          DO ifrq=1,NFREQ,IHFREQ
            FSPEC(1)=IOBSV
            FSPEC(2)=OBSV(1,IOBSV)
            FSPEC(3)=OBSV(2,IOBSV)
            FSPEC(4)=OBSV(3,IOBSV)
            if (abs(fspec(3)).lt.1.0d-15) fspec(3)=0.0d0
            if (abs(fspec(4)).lt.1.0d-15) fspec(4)=0.0d0
            FSPEC(5)=FREQ(ifrq)
            IOBFR=IOBSV+NOBSV*(ifrq-1)
            FSPEC(6)=STOKES(1,IOBFR)
            FSPEC(7)=STOKES(2,IOBFR)
            FSPEC(8)=STOKES(3,IOBFR)
            FSPEC(9)=STOKES(4,IOBFR)
            FSPEC(10)=IOBSVZ
            FSPEC(11)=IOBSVY
            FSPEC(12)=ifrq
            CALL hfm(NIDSTOK,FSPEC)
            ENDDO   !NFREQ
          ENDDO   !IOBSV

      endif !istokes

      return
      end
