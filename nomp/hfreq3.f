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

      integer ifreq,id,icycle,mfreq
      real df,flow,fhig
      real*8 fstuple(5)

      if (ipin.ne.3) return

      df=freq(2)-freq(1)
      flow=freq(1)-df/2.
      fhig=freq(nfreq)+df/2.

      if (ifreq2p.eq.1.or.freqlow.eq.freqhig) then
          DF=freqhig-freqlow
          FLOW=freqlow-DF/2.
          FHIG=freqlow+DF/2.
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
      do ifreq=1,nfreq,ihfreq
        call hfillm(id,sngl(freq(ifreq)),0.,spectot(ifreq)*1.0d-6)
      enddo   !nfreq
      call mhrout(id,icycle,' ')

      id=idfreq
      mfreq=nint((fhig-flow)/df)
      call hbook1m(id,'Photon flux through pinhole',
     &  mfreq,flow,fhig,vmx)
      do ifreq=1,nfreq,ihfreq
        call hfillm(id,sngl(freq(ifreq)),0.,wfluxt(ifreq))
      enddo   !nfreq
      call mhrout(id,icycle,' ')

      do ifreq=1,nfreq
        fstuple(1)=freq(ifreq)
        fstuple(2)=wfluxt(ifreq)
        call hfm(nidfreqp,fstuple)
      enddo   !nfreq

      if (istokes.ne.0) then

        id=icfrs0
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Mean flux-density S0 through pinhole x 1.e-6',
     &    mfreq,flow,fhig,vmx)
        do ifreq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifreq)),0.,dble(stokes(1,ifreq))*1.0d-6)
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        id=icfrs1
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Mean flux-density S1 through pinhole x 1.e-6',
     &    mfreq,flow,fhig,vmx)
        do ifreq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifreq)),0.,dble(stokes(2,ifreq))*1.0d-6)
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        call mhrout(id,icycle,' ')
        id=icfrs2
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Mean flux-density S2 through pinhole x 1.e-6',
     &    mfreq,flow,fhig,vmx)
        do ifreq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifreq)),0.,dble(stokes(3,ifreq))*1.0d-6)
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        call mhrout(id,icycle,' ')
        id=icfrs3
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Mean flux-density S3 through pinhole x 1.e-6',
     &    mfreq,flow,fhig,vmx)
        do ifreq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifreq)),0.,dble(stokes(4,ifreq))*1.0d-6)
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        id=idfrs0
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Flux S0 through pinhole',
     &    mfreq,flow,fhig,vmx)
        do ifreq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifreq)),0.,dble(wstokes(1,ifreq)))
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        id=idfrs1
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Flux S1 through pinhole',
     &    mfreq,flow,fhig,vmx)
        do ifreq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifreq)),0.,dble(wstokes(2,ifreq)))
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        call mhrout(id,icycle,' ')
        id=idfrs2
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Flux S2 through pinhole',
     &    mfreq,flow,fhig,vmx)
        do ifreq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifreq)),0.,dble(wstokes(3,ifreq)))
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        call mhrout(id,icycle,' ')
        id=idfrs3
        mfreq=nint((fhig-flow)/df)
        call hbook1m(id,'Flux S3 through pinhole',
     &    mfreq,flow,fhig,vmx)
        do ifreq=1,nfreq,ihfreq
          call hfillm(id,sngl(freq(ifreq)),0.,dble(wstokes(4,ifreq)))
        enddo   !nfreq
        call mhrout(id,icycle,' ')

        do ifreq=1,nfreq
          fstuple(1)=freq(ifreq)
          fstuple(2)=wstokes(1,ifreq)
          fstuple(3)=wstokes(2,ifreq)
          fstuple(4)=wstokes(3,ifreq)
          fstuple(5)=wstokes(4,ifreq)
          call hfm(4600,fstuple)
        enddo   !nfreq
      endif !istokes

      return
      end
