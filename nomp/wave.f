*CMZ :  4.01/02 19/04/2023  08.53.56  by  Michael Scheer
*CMZ :  4.00/17 04/11/2022  09.31.48  by  Michael Scheer
*CMZ :  4.00/16 29/09/2022  11.22.54  by  Michael Scheer
*CMZ :  4.00/13 12/11/2021  12.27.46  by  Michael Scheer
*CMZ :  4.00/07 07/06/2020  11.17.21  by  Michael Scheer
*CMZ :  3.08/01 02/04/2019  19.18.43  by  Michael Scheer
*CMZ :  3.07/01 28/03/2019  13.59.43  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  12.16.55  by  Michael Scheer
*CMZ :  3.05/01 02/05/2018  15.57.07  by  Michael Scheer
*CMZ :  3.03/04 11/07/2017  13.14.22  by  Michael Scheer
*CMZ :  3.03/02 24/11/2015  16.42.43  by  Michael Scheer
*CMZ :  3.02/03 10/11/2014  10.56.00  by  Michael Scheer
*CMZ :  3.02/00 29/08/2014  17.12.58  by  Michael Scheer
*CMZ :  3.01/06 20/06/2014  16.24.35  by  Michael Scheer
*CMZ :  3.01/05 13/06/2014  11.54.45  by  Michael Scheer
*CMZ :  3.01/00 17/06/2013  08.45.46  by  Michael Scheer
*CMZ :  3.00/01 20/03/2013  14.23.17  by  Michael Scheer
*CMZ :  2.70/06 03/01/2013  16.21.09  by  Michael Scheer
*CMZ :  2.70/00 29/11/2012  16.23.27  by  Michael Scheer
*CMZ :  2.68/02 14/06/2012  11.51.53  by  Michael Scheer
*CMZ :  2.67/04 14/05/2012  13.07.19  by  Michael Scheer
*CMZ :  2.67/02 30/04/2012  15.56.09  by  Michael Scheer
*CMZ :  2.52/16 12/11/2009  16.27.11  by  Michael Scheer
*CMZ :  2.40/00 04/03/2002  19.57.52  by  Michael Scheer
*CMZ :  2.16/08 20/10/2000  11.29.02  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  2.13/11 22/03/2000  15.25.05  by  Michael Scheer
*CMZ :  2.12/00 04/06/99  13.20.48  by  Michael Scheer
*CMZ :  2.00/00 06/01/99  11.57.36  by  Michael Scheer
*CMZ :  1.01/00 26/11/97  17.25.38  by  Michael Scheer
*CMZ : 00.01/10 29/05/96  16.48.55  by  Michael Scheer
*CMZ : 00.01/09 05/10/95  12.45.12  by  Michael Scheer
*CMZ : 00.01/04 09/12/94  11.24.52  by  Michael Scheer
*CMZ : 00.01/02 22/11/94  10.19.46  by  Michael Scheer
*CMZ : 00.01/01 23/06/94  13.26.13  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  12.02.03  by  Michael Scheer
*CMZ : 00.00/03 29/04/94  10.22.08  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.30  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WAVE
*KEEP,GPLHINT.
*KEND.

      use ompmod
      use clustermod

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.

C THIS IS THE PROCESSING ROUTINE OF THE PROGRAM WAVE
C IT CALLS THE INITIALIZATION ROUTINE GFINIT AND THE
C MODULES FOR THE INDIVIDUAL TASKS


      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,wlsopt.
      include 'wlsopt.cmn'
*KEEP,fourier.
      include 'fourier.cmn'
*KEEP,halbasy.
      include 'halbasy.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,track0.
      include 'track0.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,gseed.
      include 'gseed.cmn'
*KEEP,random.
      include 'random.cmn'
*KEEP,wvers.
      include 'wvers.cmn'
*KEND.

      INTEGER NPOL,irootmode,lun,lun99,i,lunpid,iutil_fexist

      DOUBLE PRECISION BETX0,BETY0,BETZ0,BETXF0,BETYF0,BETZF0,
     &            DTIM,BSHIFT,GAMMA

      DOUBLE PRECISION DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7,DUM8,DUM9,DUM10
      DOUBLE PRECISION DUM11,DUM12,DUM13,DUM14,DUM15,DUM16,DUM17,DUM18,DUM19,DUM20
      DOUBLE PRECISION DUM21,DUM22,DUM23,DUM24
      real dumvers
      integer idumvers

      open(newunit=lunpid,file='wave.pid')
      kpid=getpid()
      write(lunpid,*) kpid
      close(lunpid)

      chwversion=
*KEEP,wversion.
      include 'wversion.cmn'
*KEND.

      nocern=1

      do i=1,128
        if (chwversion(i:i).eq."V") then
c          print*,chwversion(i+8:i+11)
          read(chwversion(i+8:i+11),'(f4.2)')wversion
          read(chwversion(i+13:i+14),'(i2)')idumvers
          exit
        endif
      enddo

      dumvers=idumvers
      wversion=wversion+dumvers/10000.

c--- Get environment
c      call waveenvironment
c---  Check clustermod
      call waveinstances

C--- READ CONTROL-FLAGS FROM DATA-FILE AND INITIALIZE REFERENCE ORBIT

      CALL GFINIT (BETX0,BETY0,BETZ0,BETXF0,BETYF0,BETZF0,
     &             DTIM,BSHIFT,GAMMA)

C--- TERMINATE PROGRAM AFTER CALL TO UNAME

      IF (IABEND.EQ.5) goto 9999

C--- TERMINATE PROGRAM IF FOURIER COEFFICIENTS ARE WRITTEN TO FILE

      IF (IABEND.EQ.1) goto 9999

C--- TERMINATE PROGRAM ACCORDING TO WLSJUST

      IF (IABEND.EQ.2) goto 9999

C--- TERMINATE PROGRAM IF 3D POLYNOMIAL COEFFICIENTS ARE WRITTEN TO FILE

      IF (IABEND.EQ.3) THEN
c         CALL HISEND
c        RETURN
        goto 99
      ENDIF

C--- CALCULATE  A SERIES OF TRAJECTORIES AND WRITE INITIAL AND FINAL VALUES
C    ON DATA FILE

      IF (IGENFUN.GT.0) THEN
        if (iutil_fexist("wave_tranpoly.ako").ne.0) then
          print*,"*** Error in wave: File wave_tranpoly.ako already exists ***"
          print*,"*** Program WAVE aborted ***"
          print*,""
          write(lungfo,*)"*** Error in wave: File wave_tranpoly.ako already exists ***"
          write(lungfo,*)"*** Program WAVE aborted ***"
          write(lungfo,*)""
          stop
        endif
      endif

      IF (IOPTIC.NE.0) THEN
          CALL OPTI
     &            (X0,Y0,Z0,BETX0,BETY0,BETZ0,
     &            XF0,YF0,ZF0,BETXF0,BETYF0,BETZF0,
     &            DTIM,BSHIFT,GAMMA)
      ENDIF      !IOPTIC

C--- CALCULATE OPTICAL FUNCTIONS INSIDE THE DEVICE AND EMITTANCE EFFECTS

      IF (IEMIT.NE.0) CALL WBETDIS

C--- CALCULATE EMITTANCE AND POLARISATION EFFECTS OF AHW FROM ANALYTICAL
C    FORMULAS

      IF(IEMIAHW.NE.0) THEN

         CALL EMIT(   B0HALBASY,ZLHALBASY,FASYM,
     &                  DMYENERGY,RDIPOL,TAUPOL01G,BETFUN,DI2RING,DI5RING,
     &                  DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,
     &                  DUM7,DUM8,DUM9,DUM10,DUM11,DUM12,DUM13,DUM14,
     &                  DUM15,DUM16,DUM17,LUNGFO,
     &                  DISP0,DUM18,DUM19,DUM20,DUM21,DUM22,
     &                  DUM23,DUM24)


      ENDIF      !IEMIAHW

C--- UM AEQUIVALENTEN WLS ZU SUCHEN, SIEHE ALTE VERSION VON WAVE.FOR ODER
C    WLS.FOR. (SIEHE Z.B. AUF BACKUP VOM 15.7.92 ODER PROGRAM LISTINGS)
C                                                       21.7.92
C    GLEICHES GILT FUER ALTE ROUTINE DISPER.FOR

C---      SEARCH OPTIMAL PARAMETERS OF SIMPLE WLS MODEL

      IF (IWLSOPT.NE.0) THEN
                 CALL WLSOPT(
     &                1,IEMICRIT,DMYENERGY,B0MIN,B0MAX,DB0,
     &                EMICRTMX,TAUCRTMX,POLLEVMN,
     &                XLAM0MN,XLAM0MX,DXLAM0,FASYMMN,FASYMMX,DFASYM,
     &                BETFUN,RDIPOL,DBHOMF,TAUPOL01G,DI2RING,DI5RING,
     &                ZMAXMX,ZMAXMN,1,LUNGFO,IWLSADI,DISP0,DXHOM)
      ENDIF      !IWLSOPT

C--- CALCULATE COEFFICIENTS OF GENERATING FUNCTION OR MAPPING

      IF (IGENFUN.GT.0) THEN
        CALL TRANPOLY
      ELSE IF (IGENFUN.LT.0) THEN
        WRITE(6,*)'*** ERROR IN WAVE: IGENFUN.LT.0 NOT YET READY ***'
C     CALL TRANMAP ! NOCH IN ENTWICKLUNG: XF,XPF,... KANN NICHT MIT
C     EINEM SATZ VON KOEFFIZIENTEN BESCHRIBEN WERDEN, SONDERN
C     JEDE GROESSE MUSS MIT EINEM SATZ GEFITTET WERDEN. VORERST AUFGEGEBEN...
      ENDIF

C--- CALCULATE SPECTRUM

      IF(ISPEC.NE.0) CALL SPECTRUM

C--- TERMINATE PROGRAM AFTER CALL ADDAMPLI

      IF (IABEND.EQ.6) goto 9999

C--- CALCULATE POWER DENSITY ALONG BEAMLINE

      IF (IPOWER.NE.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     RESULTS FROM SR BEAMPOW:'
          WRITE(LUNGFO,*)

          CALL BEAMPOW(NPOL)

          WRITE(LUNGFO,*)
     &'     Number of poles detected:                           ',NPOL
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &'     z-Positions of beamline walls [m]:                  '
          WRITE(LUNGFO,*)
     &'                                         ',SNGL(WALL(1)),SNGL(WALL(2))
          WRITE(LUNGFO,*)
     &'     Start and end of beamline (x) [m]:                  '
          WRITE(LUNGFO,*)
     &'                                         ',SNGL(XWALLI),SNGL(XWALLE)
          WRITE(LUNGFO,*)
     &'     Minimum distance to trajectory [m]:                 ',TOTMAX(11)
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &'     x-Position of absorber:                             ',
     &      SNGL(XABSORB)
          WRITE(LUNGFO,*)
     &'     Start and end of absorber (z) [m]:                  '
          WRITE(LUNGFO,*)
     &'                                         '
     &     ,SNGL(ZABSORB(1)),SNGL(ZABSORB(2))
          WRITE(LUNGFO,*)
     &'     Minimum distance to trajectory [m]:                 ',
     &      TOTMAX(12)
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &'     Number of considered points (NPWALL):              ',NPWALL
          WRITE(LUNGFO,*)
     &'     Cutoff of magnetic field (POWBCUT):                ',
     &      SNGL(POWBCUT)
C3.5.93          WRITE(LUNGFO,*)'CORRECTION FACTOR POWCOR:',SNGL(POWCOR)
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &'     Integrated power on beamline walls [W]:            '
          WRITE(LUNGFO,*)
     &'                                         ',TOTGAM(3),TOTGAM(4)
          WRITE(LUNGFO,*)
     &'     Integrated power on absorber [W]:                  ',
     &      TOTGAM(6)

          WRITE(LUNGFO,*)
     &'     Maximum power densities on beamline walls [W/m**2]:'
          WRITE(LUNGFO,*)
     &'                                         ',TOTMAX(1),TOTMAX(2)
          WRITE(LUNGFO,*)
     &'     Maximum power density on absorber [W/m**2]:        ',TOTMAX(7)
          WRITE(LUNGFO,*)
     &'     Max. 1D power densities on beamline walls [W/m]:   '
          WRITE(LUNGFO,*)
     &'                                         ',TOTMAX(5),TOTMAX(6)
          WRITE(LUNGFO,*)
     &'     Max. 1D power density on absorber [W/m]:           ',TOTMAX(9)
          WRITE(LUNGFO,*)
     &'     Integrated photon rates on beamline walls:         '
          WRITE(LUNGFO,*)
     &'                                         ',TOTGAM(1),TOTGAM(2)
          WRITE(LUNGFO,*)
     &'     Integrated photon rate on absorber:                ',
     &      TOTGAM(5)
          WRITE(LUNGFO,*)
      ENDIF      !IPOWER

      IF (IBFORCE.NE.0) CALL BFORCE

C--- HBOOK

99      IF (IHBOOK.NE.0) then
          write(6,*)''
          write(6,*)"     Terminating histogramming"
          CALL ZEIT(6)
          write(6,*)''
          CALL HISEND
        endif

        if (bmaxgl.gt.2.5) then
          print*
          print*,"*** WARNING: LARGE B-FIELD DETECTED ***"
          print*
        endif

C--- USER OUTPUT ROUTINE

      IF (IUOUT.NE.0) CALL UOUT
      IF (IHISASCII.NE.0) CALL USERASCII

      if (iroottrees.ne.0) then
        irootmode=1
      endif

      if (iroothdf5.ne.0) then
        print*,"SR wave: iroothdf5:",iroothdf5
        print*,' '
        print*,'--- Warning: HDF5-Data only available for old versions of WAVE, supporting CERN root ---'
        print*,' '
c+self,if=linux.
c        call roottoh5
c+self,if=-linux.
c        print*,' '
c        print*,'--- Warning: HDF5-Data only available for Linux-Version of WAVE ---'
c        print*,' '
c+self. ,if=linux.
      endif


      call  util_random_get_seed(irnsize,irnseed)

      !call util_get_free_lun(lun)
      open(newunit=lun,file='wave.seeds',status='unknown')
      write(lun,*)irnsize, icode
      do i=1,irnsize
        write(lun,*)i,irnseed(i)
      enddo
      flush(lun)
      close(lun)

      RETURN

9999  iroottrees=0
      return
      END
