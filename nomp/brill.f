*CMZ :  4.01/02 12/05/2023  12.40.43  by  Michael Scheer
*CMZ :  4.00/07 16/04/2020  20.09.04  by  Michael Scheer
*CMZ :  3.05/28 18/12/2018  13.58.31  by  Michael Scheer
*CMZ :  3.05/04 27/06/2018  14.22.34  by  Michael Scheer
*CMZ :  3.05/00 25/04/2018  13.09.51  by  Michael Scheer
*CMZ :  3.02/04 21/01/2015  14.21.41  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.66/14 22/10/2010  12.10.20  by  Michael Scheer
*CMZ :  2.66/13 25/06/2010  14.49.11  by  Michael Scheer
*CMZ :  2.66/07 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.66/03 18/11/2009  10.22.04  by  Michael Scheer
*CMZ :  2.66/01 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.66/00 09/10/2009  11.37.29  by  Michael Scheer
*CMZ :  2.63/00 14/09/2009  15.19.42  by  Michael Scheer
*CMZ :  2.61/05 11/04/2007  11.58.07  by  Michael Scheer
*CMZ :  2.52/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.48/04 15/03/2004  15.22.55  by  Michael Scheer
*CMZ :  2.20/01 18/11/2000  10.19.57  by  Michael Scheer
*CMZ :  2.17/00 06/11/2000  16.13.47  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  16.27.20  by  Michael Scheer
*CMZ :  2.16/04 24/06/2000  19.53.18  by  Michael Scheer
*CMZ :  2.16/00 15/06/2000  11.56.41  by  Michael Scheer
*-- Author :    Michael Scheer   08/06/2000
      SUBROUTINE BRILL
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
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEEP,bwpolyederf90u.
      include 'bwpolyederf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,klotz.
      include 'klotz.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

C--- CONVERTION OF STOKES-VECTORS TO BRILLIANCE  (SPECTRAL BRIGHTNESS)
C--- SIGR=SQRT(2*LAMBDA*L)/2/PI*SIGRC AND SIGRP=SQRT(LAMBDA/2/L)*SIGRP
C--- THIS MAKES ONLY SENSE FOR UNDULATORS!

      DOUBLE PRECISION WAVE,SIGRP,SIGR,DEVILEN,CONV,CONVF,SIGZ,SIGY
     &                  ,SIGRF,SIGRPF,SIGZP,SIGYP,CENXEXI,DIST

      INTEGER IFREQ,IS

      IF (SIGRC.EQ.0.D0) SIGRC=1.D0
      IF (SIGRPC.EQ.0.D0) SIGRPC=1.D0

      IF (XIANF.GT.XSTART) THEN
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN BRILL: XIANF.GT.XSTART'
        WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN BRILL: XIANF.GT.XSTART'
        WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(6,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN BRILL: XIANF.GT.XSTART'
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN BRILL: XIANF.GT.XSTART'
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(LUNGFO,*)
        STOP
      ENDIF

      IF (XIEND.LT.XSTOP) THEN
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN BRILL: XIEND.LT.XSTOP'
        WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN BRILL: XIEND.LT.XSTOP'
        WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(6,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN BRILL: XIEND.LT.XSTOP'
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN BRILL: XIEND.LT.XSTOP'
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(LUNGFO,*)
        STOP
      ENDIF

      IF (NSOURCE.NE.1) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN BRILL: Number of source points not one'
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN BRILL: Number of source points not one'
        WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(6,*)
        STOP
      ENDIF

      CENXEXI=(SOURCEEO(1,1,1)+SOURCEAO(1,1,1))/2.D0
      DIST=PINCEN(1)-CENXEXI

      SIGZ=BSIGZ(ISIGSTO)
      SIGY=BSIGY(ISIGSTO)
      SIGZP=BSIGZP(ISIGSTO)
      SIGYP=BSIGYP(ISIGSTO)

      DEVILEN=XSTOP-XSTART
      if (kampli.gt.0) then
        devilen=devilen*dble(kampli)
      else IF (IAMPLI.LT.0) THEN
        DEVILEN=ABS(IAMPLI)*DEVILEN
      ENDIF

      IF (KBREC.NE.0) DEVILEN=DEVILEN-2.D0*RANGREC
      IF (KBPOLYMAG.GT.0) DEVILEN=DEVILEN-2.D0*RANGPM

      DO IFREQ=1,NFREQ

        WAVE=WTOE1/FREQ(IFREQ)*1.D-9

        ! According to Walker for sigrc=1 and sigrpc=1
        SIGR=SQRT(2.D0*WAVE*DEVILEN)/2.D0/PI1*SIGRC
        SIGRP=SQRT(WAVE/DEVILEN/2.D0)*SIGRPC

        CONV=(DIST/SIGR)**2/2.D0/PI1

        SIGRF=SQRT(SQRT(SIGR**2+SIGZ**2)*SQRT(SIGR**2+SIGY**2))
        SIGRPF=SQRT(SQRT(SIGRP**2+SIGZP**2)*SQRT(SIGRP**2+SIGYP**2))
        CONVF=(DIST/SIGRF)**2/2.D0/PI1

        BRILLC(1,IFREQ)=STOKEC(1,IFREQ)*CONV
        BRILLC(2,IFREQ)=STOKEC(2,IFREQ)*CONV
        BRILLC(3,IFREQ)=STOKEC(3,IFREQ)*CONV
        BRILLC(4,IFREQ)=STOKEC(4,IFREQ)*CONV

        if (ibunch.eq.0) then

          IF (IFOLD.NE.0) THEN

            BRILLCF(1,IFREQ)=STOKECF(1,IFREQ)*CONVF
            BRILLCF(2,IFREQ)=STOKECF(2,IFREQ)*CONVF
            BRILLCF(3,IFREQ)=STOKECF(3,IFREQ)*CONVF
            BRILLCF(4,IFREQ)=STOKECF(4,IFREQ)*CONVF

            IF (IEFOLD.NE.0) THEN
              BRILLCE(1,IFREQ)=STOKECE(1,IFREQ)*CONV
              BRILLCE(2,IFREQ)=STOKECE(2,IFREQ)*CONV
              BRILLCE(3,IFREQ)=STOKECE(3,IFREQ)*CONV
              BRILLCE(4,IFREQ)=STOKECE(4,IFREQ)*CONV
              BRILLCEF(1,IFREQ)=STOKECEF(1,IFREQ)*CONVF
              BRILLCEF(2,IFREQ)=STOKECEF(2,IFREQ)*CONVF
              BRILLCEF(3,IFREQ)=STOKECEF(3,IFREQ)*CONVF
              BRILLCEF(4,IFREQ)=STOKECEF(4,IFREQ)*CONVF
            ENDIF   !IEFOLD

          ELSE !(IFOLD.NE.0) THEN

            BRILLCF(1,IFREQ)=STOKEC(1,IFREQ)*CONVF*(SIGRP/SIGRPF)**2
            BRILLCF(2,IFREQ)=STOKEC(2,IFREQ)*CONVF*(SIGRP/SIGRPF)**2
            BRILLCF(3,IFREQ)=STOKEC(3,IFREQ)*CONVF*(SIGRP/SIGRPF)**2
            BRILLCF(4,IFREQ)=STOKEC(4,IFREQ)*CONVF*(SIGRP/SIGRPF)**2

            IF (IEFOLD.NE.0) THEN
              BRILLCE(1,IFREQ)=STOKECE(1,IFREQ)*CONV
              BRILLCE(2,IFREQ)=STOKECE(2,IFREQ)*CONV
              BRILLCE(3,IFREQ)=STOKECE(3,IFREQ)*CONV
              BRILLCE(4,IFREQ)=STOKECE(4,IFREQ)*CONV
              BRILLCEF(1,IFREQ)=STOKECE(1,IFREQ)*CONVF*(SIGRP/SIGRPF)**2
              BRILLCEF(2,IFREQ)=STOKECE(2,IFREQ)*CONVF*(SIGRP/SIGRPF)**2
              BRILLCEF(3,IFREQ)=STOKECE(3,IFREQ)*CONVF*(SIGRP/SIGRPF)**2
              BRILLCEF(4,IFREQ)=STOKECE(4,IFREQ)*CONVF*(SIGRP/SIGRPF)**2
            ENDIF   !IEFOLD

          ENDIF   !IFOLD

        endif !ibunch
      ENDDO !IFREQ

      DO IS=1,4

        IF (STOKCMX(IS,1).NE.0.D0) THEN

          BRILLCMX(IS,1)=STOKCMX(IS,1)
          WAVE=WTOE1/BRILLCMX(IS,1)*1.D-9
          SIGR=SQRT(2.D0*WAVE*DEVILEN)/2.D0/PI1*SIGRC
          SIGRP=SQRT(WAVE/DEVILEN/2.D0)*SIGRPC
          CONV=(DIST/SIGR)**2/2.D0/PI1
          SIGRF=SQRT(SQRT(SIGR**2+SIGZ**2)*SQRT(SIGR**2+SIGY**2))
          SIGRPF=SQRT(SQRT(SIGRP**2+SIGZP**2)*SQRT(SIGRP**2+SIGYP**2))
          CONVF=(DIST/SIGRF)**2/2.D0/PI1
          BRILLCMX(IS,2)=STOKCMX(IS,2)*CONV
          B3CMX(1)=G3CMX(1)
          B3CMX(2)=G3CMX(2)*CONV

c Estimated flux
          BRILLCMX(IS+4,2)=STOKCMX(IS,2)*dist**2*2.0d0*pi1*sigrp**2

          if (ibunch.eq.0) then
            IF (IEFOLD.NE.0.and.STOKCMXE(IS,1).ne.0.0d0) THEN
              BRILLCMXE(IS,1)=STOKCMXE(IS,1)
              WAVE=WTOE1/BRILLCMXE(IS,1)*1.D-9
              SIGR=SQRT(2.D0*WAVE*DEVILEN)/2.D0/PI1*SIGRC
              SIGRP=SQRT(WAVE/DEVILEN/2.D0)*SIGRPC
              CONV=(DIST/SIGR)**2/2.D0/PI1
              SIGRF=SQRT(SQRT(SIGR**2+SIGZ**2)*SQRT(SIGR**2+SIGY**2))
              SIGRPF=SQRT(SQRT(SIGRP**2+SIGZP**2)*SQRT(SIGRP**2+SIGYP**2))
              CONVF=(DIST/SIGRF)**2/2.D0/PI1
              BRILLCMXE(IS,2)=STOKCMXE(IS,2)*CONV
c Estimated flux
              BRILLCMXE(IS+4,2)=STOKCMXE(IS,2)*dist**2*2.0d0*pi1*sigrp**2
              B3CMXE(1)=G3CMXE(1)
              B3CMXE(2)=G3CMXE(2)*CONV
            ENDIF  !IEFOLD

            IF (IFOLD.NE.0.and.STOKCMXF(IS,1).ne.0.0d0) THEN
              BRILLCMXF(IS,1)=STOKCMXF(IS,1)
              WAVE=WTOE1/BRILLCMXF(IS,1)*1.D-9
              SIGR=SQRT(2.D0*WAVE*DEVILEN)/2.D0/PI1*SIGRC
              SIGRP=SQRT(WAVE/DEVILEN/2.D0)*SIGRPC
              CONV=(DIST/SIGR)**2/2.D0/PI1
              SIGRF=SQRT(SQRT(SIGR**2+SIGZ**2)*SQRT(SIGR**2+SIGY**2))
              SIGRPF=SQRT(SQRT(SIGRP**2+SIGZP**2)*SQRT(SIGRP**2+SIGYP**2))
              CONVF=(DIST/SIGRF)**2/2.D0/PI1
              BRILLCMXF(IS,2)=STOKCMXF(IS,2)*CONVF
c Estimated flux
              BRILLCMXF(IS+4,2)=STOKCMXF(IS,2)*dist**2*2.0d0*pi1*sigrpf**2
              B3CMXF(1)=G3CMXF(1)
              B3CMXF(2)=G3CMXF(2)*CONVF
              IF (IEFOLD.NE.0.and.STOKCMXEF(IS,1).ne.0.0d0) THEN
                BRILLCMXEF(IS,1)=STOKCMXEF(IS,1)
                WAVE=WTOE1/BRILLCMXEF(IS,1)*1.D-9
                SIGR=SQRT(2.D0*WAVE*DEVILEN)/2.D0/PI1*SIGRC
                SIGRP=SQRT(WAVE/DEVILEN/2.D0)*SIGRPC
                CONV=(DIST/SIGR)**2/2.D0/PI1
                SIGRF=SQRT(SQRT(SIGR**2+SIGZ**2)*SQRT(SIGR**2+SIGY**2))
                SIGRPF=SQRT(SQRT(SIGRP**2+SIGZP**2)*SQRT(SIGRP**2+SIGYP**2))
                CONVF=(DIST/SIGRF)**2/2.D0/PI1
                BRILLCMXEF(IS,2)=STOKCMXEF(IS,2)*CONVF
c Estimated flux
                BRILLCMXEF(IS+4,2)=STOKCMXEF(IS,2)*dist**2*2.0d0*pi1*sigrp**2
                B3CMXEF(1)=G3CMXEF(1)
                B3CMXEF(2)=G3CMXEF(2)*CONVF
              ENDIF   !IEFOLD
            ELSE IF (STOKCMX(IS,1).ne.0.0d0) THEN   !IFOLD
              BRILLCMXF(IS,1)=STOKCMX(IS,1)
              WAVE=WTOE1/BRILLCMX(IS,1)*1.D-9
              SIGR=SQRT(2.D0*WAVE*DEVILEN)/2.D0/PI1*SIGRC
              SIGRP=SQRT(WAVE/DEVILEN/2.D0)*SIGRPC
              CONV=(DIST/SIGR)**2/2.D0/PI1
              SIGRF=SQRT(SQRT(SIGR**2+SIGZ**2)*SQRT(SIGR**2+SIGY**2))
              SIGRPF=SQRT(SQRT(SIGRP**2+SIGZP**2)*SQRT(SIGRP**2+SIGYP**2))
              CONVF=(DIST/SIGRF)**2/2.D0/PI1
              STOKCMXF(IS,2)=STOKCMX(IS,2)*(SIGRP/SIGRPF)**2
              BRILLCMXF(IS,2)=STOKCMX(IS,2)*CONVF*(SIGRP/SIGRPF)**2
c Estimated flux
              ! flux and folded flux are the same, therefore sigrp!!
              BRILLCMXF(IS+4,2)=STOKCMX(IS,2)*dist**2*2.0d0*pi1*sigrp**2
              B3CMXF(1)=G3CMX(1)
              B3CMXF(2)=G3CMX(2)*CONVF*(SIGRP/SIGRPF)**2
              STOKCMXEF(IS,2)=STOKCMXE(IS,2)*(SIGRP/SIGRPF)**2
              IF (IEFOLD.NE.0.and.STOKCMXE(IS,1).ne.0.0d0) THEN
                BRILLCMXEF(IS,1)=STOKCMXE(IS,1)
                WAVE=WTOE1/BRILLCMXE(IS,1)*1.D-9
                SIGR=SQRT(2.D0*WAVE*DEVILEN)/2.D0/PI1*SIGRC
                SIGRP=SQRT(WAVE/DEVILEN/2.D0)*SIGRPC
                CONV=(DIST/SIGR)**2/2.D0/PI1
                SIGRF=SQRT(SQRT(SIGR**2+SIGZ**2)*SQRT(SIGR**2+SIGY**2))
                SIGRPF=SQRT(SQRT(SIGRP**2+SIGZP**2)*SQRT(SIGRP**2+SIGYP**2))
                CONVF=(DIST/SIGRF)**2/2.D0/PI1
                BRILLCMXEF(IS,2)=STOKCMXE(IS,2)*CONVF*(SIGRP/SIGRPF)**2
c Estimated flux
                ! flux and folded flux are the same, therefore sigrp!!
                BRILLCMXEF(IS+4,2)=STOKCMXE(IS,2)*dist**2*2.0d0*pi1*sigrp**2
                B3CMXEF(1)=G3CMXE(1)
                B3CMXEF(2)=G3CMXE(2)*CONVF*(SIGRP/SIGRPF)**2
              ENDIF   !IEFOLD
            ENDIF  !IFOLD
          ENDIF !STOKCMX
        endif !ibunch
      ENDDO  !IS

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

      write(lungfo,*)'     Brilliance calculations:'
      write(lungfo,*)'     (only meaningful for undulators)'
      write(lungfo,*)
      write(lungfo,*)
     &  '     sigr=sqrt(2*lambda*L)/2/pi*sigrc and '
      write(lungfo,*)
     &  '     sigrp=sqrt(lambda/2/L)*sigrpc'
      write(lungfo,*)
      write(lungfo,*)'      sigrc,sigrpc:',sigrc,sigrpc
      write(lungfo,*)
      write(lungfo,*)
     &  '     make sure that the device length agrees to the length'
      write(lungfo,*)
     &  '     of the trajectory'
      write(lungfo,*)'     (for rec-fields rangrec and rangpm are substracted by program)'

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     Observation point:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     '
     &  ,SNGL(OBSV(1,ICBRILL)),SNGL(OBSV(2,ICBRILL)),SNGL(OBSV(3,ICBRILL))
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     considered device length [m]:',SNGL(DEVILEN)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     Brilliance [photons/sec/BW/mm**2/mrad**2]'
      WRITE(LUNGFO,*)

      DO IFREQ=1,NFREQ
        WRITE(LUNGFO,'(5x,5G15.5)')SNGL(FREQ(IFREQ))
     &    ,(BRILLC(IS,IFREQ)*(1.D-6)**2,IS=1,4)
      ENDDO !IFREQ

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     estimated brilliance maxima, flux and flux-density of S0, S1, S2, S3:'
      WRITE(LUNGFO,*)
      DO IS=1,4
        WRITE(LUNGFO,'(5X,4G15.5)')BRILLCMX(IS,1)
     &    ,BRILLCMX(IS,2)*(1.D-6)**2
     &    ,BRILLCMX(IS+4,2)
     &    ,STOKCMX(IS,2)*DIST**2*(1.D-3)**2
      ENDDO

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     estimated brilliance maximum of S3*S3/S0:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,'(5X,2G15.5)')B3CMX(1),B3CMX(2)*(1.D-6)**2

      if (ibunch.eq.0) then

        IF (IFOLD.EQ.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '     ---> flag IFOLD is not set, emittance effects are estimated'
          WRITE(LUNGFO,*)
     &      '     from horizontal beam size [mm] and divergence [mrad]: '
          WRITE(LUNGFO,*)
     &      '     ',SNGL(BSIGZ(ISIGSTO)),SNGL(BSIGZP(ISIGSTO))
          WRITE(LUNGFO,*)
     &      '     vertical beam size [mm] and divergence [mrad]: '
          WRITE(LUNGFO,*)
     &      '     ',SNGL(BSIGY(ISIGSTO)),SNGL(BSIGYP(ISIGSTO))
          WRITE(LUNGFO,*)
        ENDIF   !IFOLD

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     Brilliance [photons/sec/BW/mm**2/mrad**2] with emittance'
        WRITE(LUNGFO,*)

        DO IFREQ=1,NFREQ
          WRITE(LUNGFO,'(5x,5G15.5)')SNGL(FREQ(IFREQ))
     &      ,(BRILLCF(IS,IFREQ)*(1.D-6)**2,IS=1,4)
        ENDDO   !IFREQ

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     estimated brilliance maxima, flux and flux-density of S0, S1, S2, S3'
        WRITE(LUNGFO,*)
     &    '     with emittance:'
        WRITE(LUNGFO,*)
        DO IS=1,4
          WRITE(LUNGFO,'(5X,4G15.5)')BRILLCMXF(IS,1)
     &      ,BRILLCMXF(IS,2)*(1.D-6)**2
     &      ,BRILLCMXF(IS+4,2)
     &      ,STOKCMXF(IS,2)*DIST**2*(1.D-3)**2
        ENDDO

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     estimated brilliance maximum of S3*S3/S0 with emittance:'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,'(5X,2G15.5)')B3CMXF(1),B3CMXF(2)*(1.D-6)**2

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     corresponding sigmas of source size [mm]'
        WRITE(LUNGFO,*)
     &    '     and divergence [mrad] without and with emittance'
        WRITE(LUNGFO,*)
     &    '     (including corr. factors sigrc,sigrpc):'
        IF (BRILLCMXF(1,1).NE.0.0D0) THEN
          WAVE=WTOE1/BRILLCMXF(1,1)*1.D-9
        ELSE
          WAVE=0.0D0
        ENDIF
        SIGR=SQRT(2.D0*WAVE*DEVILEN)/2.D0/PI1*SIGRC
        SIGRP=SQRT(WAVE/DEVILEN/2.D0)*SIGRPC
        SIGRF=SQRT(SQRT(SIGR**2+SIGZ**2)*SQRT(SIGR**2+SIGY**2))
        SIGRPF=SQRT(SQRT(SIGRP**2+SIGZP**2)*SQRT(SIGRP**2+SIGYP**2))
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     L/mm, lambda/mm: ',SNGL(DEVILEN*1000.0d0),SNGL(WAVE*1000.0d0)
        WRITE(LUNGFO,*)'     without emit.: ',SNGL(SIGR*1000.D0),SNGL(SIGRP*1000.D0)
        WRITE(LUNGFO,*)'        with emit.: ',SNGL(SIGRF*1000.D0),SNGL(SIGRPF*1000.D0)
        WRITE(LUNGFO,*)'        (sigrf**2 = sigrxf*sigryf)'
        WRITE(LUNGFO,*)

        IF (IEFOLD.NE.0) THEN

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '     Brilliance [photons/sec/BW/mm**2/mrad**2] (e-folded)'
          WRITE(LUNGFO,*)

          DO IFREQ=1,NFREQ
            WRITE(LUNGFO,'(5x,5G15.5)')SNGL(FREQ(IFREQ))
     &        ,(BRILLCE(IS,IFREQ)*(1.D-6)**2,IS=1,4)
          ENDDO !IFREQ

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     estimated brilliance maxima, flux and flux-density of S0, S1, S2, S3'
          WRITE(LUNGFO,*)'     with energy spread:'
          WRITE(LUNGFO,*)
          DO IS=1,4
            WRITE(LUNGFO,'(5X,4G15.5)')BRILLCMXE(IS,1)
     &        ,BRILLCMXE(IS,2)*(1.D-6)**2
     &        ,BRILLCMXE(IS+4,2)
     &        ,STOKCMXE(IS,2)*DIST**2*(1.D-3)**2
          ENDDO

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     estimated brilliance maximum of S3*S3/S0'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     with energy spread:'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,'(5X,2G15.5)')B3CMXE(1),B3CMXE(2)*(1.D-6)**2

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '     Brilliance [photons/sec/BW/mm**2/mrad**2] with emittance'
          WRITE(LUNGFO,*)
     &      '     (e-folded)'
          WRITE(LUNGFO,*)

          DO IFREQ=1,NFREQ
            WRITE(LUNGFO,'(5x,5G15.5)')SNGL(FREQ(IFREQ))
     &        ,(BRILLCEF(IS,IFREQ)*(1.D-6)**2,IS=1,4)
          ENDDO !IFREQ

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '     estimated brilliance maxima, flux and flux-density of S0, S1, S2, S3'
          WRITE(LUNGFO,*)
     &      '     with emittance and energy spread:'
          WRITE(LUNGFO,*)
          DO IS=1,4
            WRITE(LUNGFO,'(5X,4G15.5)')BRILLCMXEF(IS,1)
     &        ,BRILLCMXEF(IS,2)*(1.D-6)**2
     &        ,BRILLCMXEF(IS+4,2)
     &        ,STOKCMXEF(IS,2)*DIST**2*(1.D-3)**2
          ENDDO

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '     estimated brilliance maximum of S3*S3/S0'
          WRITE(LUNGFO,*)
     &      '     with emittance and energy spread:'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,'(5X,2G15.5)')B3CMXEF(1),B3CMXEF(2)*(1.D-6)**2

        ENDIF   !EFOLD

      endif !ibunch

      IF(IWFILBRILL.NE.0.AND.ISTOKES.NE.0) THEN

        OPEN(UNIT=LUNC,FILE=FILEBRILL,STATUS='NEW')

        WRITE(LUNC,*)ICODE,' ',CODE
        WRITE(LUNC,*)NFREQ

        DO IFREQ=1,NFREQ
          IF (IUNIT.EQ.0)
     &      WRITE(LUNC,*)SNGL(FREQ(IFREQ))
     &      ,((BRILLC(IS,IFREQ)),IS=1,4)
        ENDDO

        CLOSE(LUNC)
      ENDIF !IWFILB

      if (ibunch.eq.0) then

        IF(IWFILBRILLF.NE.0.AND.ISTOKES.NE.0) THEN

          OPEN(UNIT=LUNCF,FILE=FILEBRILLF,STATUS='NEW')

          WRITE(LUNCF,*)ICODE,' ',CODE
          WRITE(LUNCF,*)NFREQ

          DO IFREQ=1,NFREQ
            IF (IUNIT.EQ.0)
     &        WRITE(LUNCF,*)SNGL(FREQ(IFREQ))
     &        ,((BRILLCF(IS,IFREQ)),IS=1,4)
          ENDDO

          CLOSE(LUNCF)
        ENDIF !IWFILBF

        IF(IEFOLD.NE.0.AND.IWFILBRILLE.NE.0.AND.ISTOKES.NE.0) THEN

          OPEN(UNIT=LUNCE,FILE=FILEBRILLE,STATUS='NEW')

          WRITE(LUNCE,*)ICODE,' ',CODE
          WRITE(LUNCE,*)NFREQ

          DO IFREQ=1,NFREQ
            IF (IUNIT.EQ.0)
     &        WRITE(LUNCE,*)SNGL(FREQ(IFREQ))
     &        ,((BRILLCE(IS,IFREQ)),IS=1,4)
          ENDDO
          CLOSE(LUNCE)
        ENDIF !IWFILBE

        IF(IEFOLD.NE.0.AND.IWFILBRILLEF.NE.0.AND.ISTOKES.NE.0) THEN

          OPEN(UNIT=LUNCEF,FILE=FILEBRILLEF,STATUS='NEW')

          WRITE(LUNCEF,*)ICODE,' ',CODE
          WRITE(LUNCEF,*)NFREQ

          DO IFREQ=1,NFREQ
            IF (IUNIT.EQ.0)
     &        WRITE(LUNCEF,*)SNGL(FREQ(IFREQ))
     &        ,((BRILLCEF(IS,IFREQ)),IS=1,4)
          ENDDO

          CLOSE(LUNCEF)

        ENDIF !IWFILBEF

      endif !ibunch

      OPEN(UNIT=LUNCEF,FILE='brill_brilliance.dat',
     &  STATUS='unknown')
      do is=1,4
        WRITE(LUNcef,'(2G15.5)')BRILLCMX(IS,1)
     &    ,BRILLCMX(IS,2)*(1.D-6)**2
      enddo
      close(luncef)

      OPEN(UNIT=LUNCEF,FILE='brill_flux-density.dat',
     &  STATUS='unknown')
      do is=1,4
        WRITE(LUNcef,'(2G15.5)')BRILLCMX(IS,1)
     &    ,STOKCMX(IS,2)*DIST**2*(1.D-3)**2
      enddo
      close(luncef)

      OPEN(UNIT=LUNCEF,FILE='brill_flux.dat',
     &  STATUS='unknown')
      do is=1,4
        WRITE(LUNcef,'(2G15.5)')BRILLCMX(IS,1)
     &    ,BRILLCMX(IS+4,2)
      enddo
      close(luncef)

      if (ibunch.eq.0) then

        OPEN(UNIT=LUNCEF,FILE='brill_brilliance_e.dat',
     &    STATUS='unknown')
        do is=1,4
          WRITE(LUNcef,'(2G15.5)')BRILLCMXe(IS,1)
     &      ,BRILLCMXe(IS,2)*(1.D-6)**2
        enddo
        close(luncef)

        OPEN(UNIT=LUNCEF,FILE='brill_flux-density_e.dat',
     &    STATUS='unknown')
        do is=1,4
          WRITE(LUNcef,'(2G15.5)')BRILLCMXE(IS,1)
     &      ,STOKCMXE(IS,2)*DIST**2*(1.D-3)**2
        enddo
        close(luncef)

        OPEN(UNIT=LUNCEF,FILE='brill_flux_e.dat',
     &    STATUS='unknown')
        do is=1,4
          WRITE(LUNcef,'(2G15.5)')BRILLCMXe(IS,1)
     &      ,BRILLCMXe(IS+4,2)
        enddo
        close(luncef)

        OPEN(UNIT=LUNCEF,FILE='brill_brilliance_f.dat',
     &    STATUS='unknown')
        do is=1,4
          WRITE(LUNcef,'(2G15.5)')BRILLCMXf(IS,1)
     &      ,BRILLCMXf(IS,2)*(1.D-6)**2
        enddo
        close(luncef)

        OPEN(UNIT=LUNCEF,FILE='brill_flux-density_f.dat',
     &    STATUS='unknown')
        do is=1,4
          WRITE(LUNcef,'(2G15.5)')BRILLCMXf(IS,1)
     &      ,STOKCMXf(IS,2)*DIST**2*(1.D-3)**2
        enddo
        close(luncef)

        OPEN(UNIT=LUNCEF,FILE='brill_flux_f.dat',
     &    STATUS='unknown')
        do is=1,4
          WRITE(LUNcef,'(2G15.5)')BRILLCMXf(IS,1)
     &      ,BRILLCMXf(IS+4,2)
        enddo
        close(luncef)

        OPEN(UNIT=LUNCEF,FILE='brill_brilliance_ef.dat',
     &    STATUS='unknown')
        do is=1,4
          WRITE(LUNcef,'(2G15.5)')BRILLCMXef(IS,1)
     &      ,BRILLCMXef(IS,2)*(1.D-6)**2
        enddo
        close(luncef)

        OPEN(UNIT=LUNCEF,FILE='brill_flux-density_ef.dat',
     &    STATUS='unknown')
        do is=1,4
          WRITE(LUNcef,'(2G15.5)')BRILLCMXef(IS,1)
     &      ,STOKCMXef(IS,2)*DIST**2*(1.D-3)**2
        enddo
        close(luncef)

        OPEN(UNIT=LUNCEF,FILE='brill_flux_ef.dat',
     &   STATUS='unknown')
        do is=1,4
          WRITE(LUNcef,'(2G15.5)')BRILLCMXef(IS,1)
     &      ,BRILLCMXef(IS+4,2)
        enddo
        close(luncef)

      endif !ibunch

      RETURN
      END
