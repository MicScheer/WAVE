*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.66/07 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.66/06 18/11/2009  13.56.29  by  Michael Scheer
*CMZ :  2.66/03 28/10/2009  13.07.43  by  Michael Scheer
*CMZ :  2.66/01 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.65/02 29/09/2009  09.46.33  by  Michael Scheer
*CMZ :  2.65/01 21/09/2009  14.35.04  by  Michael Scheer
*CMZ :  2.65/00 18/09/2009  07.26.13  by  Michael Scheer
*CMZ :  2.64/07 16/09/2009  12.42.25  by  Michael Scheer
*CMZ :  2.64/06 15/09/2009  15.03.21  by  Michael Scheer
*CMZ :  2.64/05 14/09/2009  09.10.19  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine pininr
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
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

C--- INITIALIZE CYLINDRICAL GRID OF OBERSERVATION POINTS OF PINHOLE

      implicit none

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      double precision r,phi
      integer iob,iphi,ir

      if (mpinr.eq.1) then
        mpinr=2
        write(6,*)
        write(6,*)
     &    '*** Warning in PININR: MPINR set form 1 to minimum value of 2'
        write(6,*)
        write(lungfo,*)
        write(lungfo,*)
     &    '*** Warning in PININR: MPINR set form 1 to minimum value of 2'
        write(lungfo,*)
      endif

      if (pinrad.eq.0.0d0) pinrad=pinr

      if (pinrad.eq.9999.0d0) then
        pinrad=-1.0d0
        do iob=1,nobsv
          r=sqrt(obsv(2,iob)**2+obsv(3,iob)**2)
          if (r.gt.pinrad) pinrad=r
        enddo
      endif

      if (abs(mpinr).eq.9999) then
        mpinr=(sqrt(2.)*max(mpinz,mpiny)/2)*2+1
      endif

      if (obsvdr.eq.9999.0d0) then
        obsvdr=min(obsvdz,obsvdy)
      endif

      if (obsvdr.ne.0.0d0) then
        mpinr=nint(pinrad/obsvdr)+2
        pinrad=obsvdr*(mpinr-1)
      else if (mpinr.gt.1) then
        obsvdr=pinrad/(mpinr-1)
      else
        stop '*** Error in PININR: Bad MPINR or OBSVDR!'
      endif !(obsvdz.ne.0.d0)

      if (pinrad.lt.sqrt((pinw/2.0d0)**2+(pinh/2.0d0)**2)) then
        print*,
     &    '*** Error in PININR: PINRAD.LT.SQRT((PINW/2.0D0)**2+(PINH/2.0d0)**2)'
        print*,'*** Please check PINRAD, OBSVDR etc.'
        stop '*** Program WAVE aborted ***'
      endif

C--- DATA OF PINHOLE ARE TAKEN FORM NAMELIST

      obsvdphi=obsvdphi/360.0d0*twopi1

      if (iquadphi.ne.0.and.(istokes.ne.0.or.iphase.ne.0)) then
c        iquadphi=0
        write(lungfo,*)' '
        write(lungfo,*)'*** Warning in PININR: IQUADPHI and ISTOKES/IPHASE are incompatible!'
        write(lungfo,*)'*** In general, only S0 will be correct'
c        write(lungfo,*)'*** Warning in PININR: IQUADPHI set to zero'
        write(6,*)' '
        write(6,*)'*** Warning in PININR: IQUADPHI and ISTOKES/IPHASE are incompatible!'
        write(6,*)'*** In general, only S0 will be correct'
c        write(6,*)'*** Warning in PININR: IQUADPHI set to zero'
      endif

      if (iquadphi.eq.0) then
        if (obsvdphi.ne.0.0d0) then
          mpinphi=nint(twopi1/obsvdphi)
        else
          obsvdphi=twopi1/mpinphi
        endif
      else
        if (obsvdphi.ne.0.0d0) then
          mpinphi=nint(pi1/2.0d0/obsvdphi)+1
        else if (mpinphi.gt.1) then
          obsvdphi=pi1/2.0d0/(mpinphi-1)
        else
          obsvdphi=twopi1
        endif
      endif

      IF (IF1DIM.NE.0.AND.MPINR.NE.1) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN PININR ***'
        WRITE(LUNGFO,*)'FLAG IF1DIM SET BUT'
        WRITE(LUNGFO,*)'MPINZ NOT EQUAL ONE, ADJUSTED'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN PININR ***'
        WRITE(6,*)'FLAG IF1DIM SET BUT'
        WRITE(6,*)'MPINZ NOT EQUAL ONE, ADJUSTED'
        WRITE(6,*)
        MPINR=1
      ENDIF

      if (iusem.ne.0) then
        stop '*** error in pininr: iusem not allowed here! '
      endif !iusem

      nobsvr=mpinr
      nobsvphi=mpinphi
      nobsvrphi=nobsvr*nobsvphi

      if (nobsvrphi.le.0) then
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN PININR ***'
        WRITE(LUNGFO,*)'Number of observation points not positive!'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN PININR ***'
        WRITE(6,*)'Number of observation points not positive!'
        WRITE(6,*)
        STOP '*** PROGRAM WAVE ABORTED ***'
      ENDIF

      if (nobsvrphi.le.0) then
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN PININR ***'
        WRITE(LUNGFO,*)'Number of observation points not positive!'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN PININR ***'
        WRITE(6,*)'Number of observation points not positive!'
        WRITE(6,*)
        STOP '*** PROGRAM WAVE ABORTED ***'
      ENDIF

c- increase pinhole; size of rectangular pinhole is already calculated in PININ

      mobsvr=nobsvr  !store values
      mobsvphi=nobsvphi
      mobsvrphi=nobsvrphi

      if (obsvdz*(nobsvz-1)/2.0d0.gt.(obsvdr*(nobsvr-1))) then
        nobsvr=obsvdz*(nobsvz-1)/2.0d0/obsvdr+1
      endif

      if (obsvdy*(nobsvy-1)/2.0d0.gt.(obsvdr*(nobsvr-1))) then
        nobsvr=obsvdy*(nobsvy-1)/2.0d0/obsvdr+1
      endif

      nobsvrphi=nobsvr*nobsvphi

      if (iobsvrphi_a.ne.nobsv) then
        if (iobsvrphi_a.ne.0) deallocate(obsvrphi)
        allocate(obsvrphi(3,nobsvrphi))
        iobsvrphi_a=nobsvrphi
      endif !(iobsv_a.lt.nobsv)

      if (iobsvr_a.ne.nobsvr) then
        if (iobsvr_a.ne.0) deallocate(obsvr)
        allocate(obsvr(nobsvr))
        iobsvr_a=nobsvr
      endif !(iobsvy_a.lt.nobsvy)

      if (iobsvphi_a.ne.nobsvphi) then
        if (iobsvphi_a.ne.0) deallocate(obsvphi)
        allocate(obsvphi(nobsvphi))
        iobsvphi_a=nobsvphi
      endif !(iobsv_a.lt.nobsv)

      iob=0
      do iphi=1,nobsvphi
        do ir=1,nobsvr
          iob=iob+1
          obsvrphi(1,iob)=pincen(1)
          r=obsvdr*(ir-1)
          phi=obsvdphi*(iphi-1)
          obsvrphi(2,iob)=r
          obsvrphi(3,iob)=phi
        enddo
      enddo

      do ir=1,nobsvr
        obsvr(ir)=obsvrphi(2,ir)
      enddo

      do iphi=1,nobsvphi
        iob=(iphi-1)*nobsvr+1
        obsvphi(iphi)=obsvrphi(3,iob)
      enddo

      write(lungfo,*)
      write(lungfo,*)'      Subroutine PININR:'
      write(lungfo,*)
      write(lungfo,*)'      NOBSVR, PINRAD, OBSVDR:',
     &  MPINR, SNGL(PINR), SNGL(OBSVDR)
      write(lungfo,*)'      IQUADPHI, NOBSVPHI, OBSVDPHI:',
     &  IQUADPHI,NOBSVPHI,SNGL(OBSVDPHI)
      write(lungfo,*)

      return
      end
