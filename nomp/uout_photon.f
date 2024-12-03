*CMZ :          17/10/2023  13.35.53  by  Michael Scheer
*CMZ :  2.15/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ : 00.02/05 19/03/97  14.10.43  by  Michael Scheer
*CMZ : 00.02/04 26/02/97  10.23.26  by  Michael Scheer
*-- Author :    Michael Scheer   25/02/97

      SUBROUTINE UOUT_photon

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.

C     INTERFACE FUER PHOTON
C     USER(1) MUSS KRITISCHE ENERGIE ENTHALTEN

C     EIGENTLICH IST DIE KRITISCHE ENERGIE UND GAMMA UNWICHTIG,
C     ABER VORERST FUER CHECKS MITAUSGEBEN.
C     PINCEN(2) MUSS NULL SEIN (FAKTOR 2)
C
C     AUSGEBEN WIRD dFLUX/dTHETA[mrad]

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,observ.
      include 'observ.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEND.


      INTEGER IFREQ
      DOUBLE PRECISION EGAMMA,PHOTONS

      WRITE(6,*)
      WRITE(6,*)'SUBROUTINE UOUT_PHOTON:'
      WRITE(6,*)'======================='
      WRITE(6,*)
      WRITE(6,*)'Writing file WAVE.PHOTON for program PHOTON'
      WRITE(6,*)'Ecrit:',ecphoton
      WRITE(6,*)


      IF (IPIN.NE.1) STOP '*** ERROR IN UOUT_PHOTON: IPIN.NE.1'
      IF (IF1DIM.NE.1) STOP '*** ERROR IN UOUT_PHOTON: IF1DIM.NE.1'
      IF (PINCEN(2).NE.0.0)STOP
     &  '*** ERROR IN UOUT_PHOTON: PINCEN(2).NE.0.0'

      OPEN(UNIT=99,FILE='WAVE.PHOTON',STATUS='NEW')

      WRITE(99,*)'WAVE.PHOTON'
      WRITE(99,*)ICODE,' ',CODE
      WRITE(99,*)ECPHOTON,DMYGAMMA
      WRITE(99,*)PINW,PINH
      WRITE(99,*)OBSVDZ,OBSVDY
      WRITE(99,*)PINCEN
      WRITE(99,*)NFREQ

      DO IFREQ=1,NFREQ
        EGAMMA=FREQ(IFREQ)/ECPHOTON
        PHOTONS=WFLUXT(IFREQ)
     &    /1.327D13 ! 1.327e13 * Ebeam**2  * Icurr * H2(y) yields flux-dens.
     &    /BANWID*0.001D0
     &    *(PINCEN(1)/PINW/1000.0D0)
     &    /(DMYCUR*1000.0D0)
     &    /DMYENERGY**2
     &    /2.D0
        WRITE(99,*)EGAMMA,PHOTONS
      ENDDO

      CLOSE(99)

      RETURN
      END
