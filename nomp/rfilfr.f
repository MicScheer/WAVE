*CMZ :  4.00/15 07/04/2022  07.14.03  by  Michael Scheer
*CMZ :  4.00/04 10/05/2019  16.57.06  by  Michael Scheer
*CMZ :  3.03/04 03/08/2017  14.10.45  by  Michael Scheer
*CMZ :  3.03/02 03/12/2015  17.38.18  by  Michael Scheer
*CMZ :  3.02/07 03/12/2015  17.37.47  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.69/00 26/10/2012  13.26.40  by  Michael Scheer
*CMZ :  2.67/02 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.66/12 22/05/2010  16.51.50  by  Michael Scheer
*CMZ :  2.66/09 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.52/01 18/11/2009  10.22.04  by  Michael Scheer
*CMZ :  2.52/00 24/06/2004  17.38.38  by  Michael Scheer
*CMZ :  2.20/10 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.20/01 21/11/2000  19.24.54  by  Michael Scheer
*CMZ :  2.16/08 25/10/2000  12.21.53  by  Michael Scheer
*CMZ :  2.16/07 14/09/2000  17.12.07  by  Michael Scheer
*CMZ :  2.16/06 30/08/2000  13.39.46  by  Michael Scheer
*CMZ :  2.16/05 02/08/2000  13.39.56  by  Michael Scheer
*CMZ :  2.16/04 21/07/2000  10.14.36  by  Michael Scheer
*CMZ :  2.15/00 03/05/2000  16.31.03  by  Michael Scheer
*CMZ : 00.02/00 19/11/96  14.56.43  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.14.45  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.53.37  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.43  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE RFILFR
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

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

C--- SUBROUTINE READS FREQUENCES FOR WHICH PHOTON FLUX IS CALCULATED

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,klotz.
      include 'klotz.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER IFR,I,J

C260194 {
      DOUBLE PRECISION DFREQ,DUM
C260194 }
C24.394 {
      INTEGER INCREAS,NFREQO,IFREQ
      DOUBLE PRECISION DEXTEND,PHIDEFL,B0DIP,DEFLEC,DLENG,HARM,DNPER,HARMWIDTH
      DOUBLE PRECISION WAVE,df
C24.394 }

      double precision frscaleo,freql,freqh

      frscaleo=frscale
      frscale=abs(frscale)

      if (pinwsc.le.0.0d0) pinwsc=1.0d0
      if (pinhsc.le.0.0d0) pinhsc=1.0d0

      IF (FRSCALE.EQ.0.D0) FRSCALE=1.D0

      PHIDEFL=ANGRMS*SQRT(2.)
      B0DIP=BRMS*SQRT(2.)
      DEFLEC=PHIDEFL*DMYGAMMA

      IF (B0DIP.NE.0.D0) THEN
        DLENG=DEFLEC/93.4/B0DIP
      ELSE
        DLENG=0.D0
      ENDIF

      IF (DLENG.EQ.0.0) DLENG=(XMX-XMN)

      HARM=0.95D3*DMYENERGY**2/(1.D0+DEFLEC**2/2.D0)/DLENG/100.d0

      IF (KBREC.EQ.0) THEN
        DNPER=(XMX-XMN)/DLENG
      ELSE  !KBREC
        DNPER=(XMX-XMN-2.D0*RANGREC)/DLENG
      ENDIF !KBREC

      IF (DNPER.LT.1.D0) THEN
        DNPER=1.D0
        HARMWIDTH=HARM*0.5
      ELSE
        HARMWIDTH=HARM/DNPER
      ENDIF

      WAVE=WTOE1/HARM*1.D-9/FRSCALE

      IF (IUNDULATOR.NE.0) THEN
        IF (FREQLOW.EQ.-9999.) THEN
          FREQLOW=HARM-HARMWIDTH
        ENDIF   !FREQLOW
        IF (FREQHIG.EQ.-9999.) THEN
          FREQHIG=HARM+HARMWIDTH
        ENDIF   !FREQLOW
        IF (PINW.EQ.9999..AND.IPIN.NE.0) THEN
          IF (IAMPLI.GT.0) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
     &        '*** ERROR IN RFILRF: PINW=9999. AND IAMPLI GREATER THAN ZERO***'
            WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
            WRITE(6,*)
            WRITE(6,*)
     &        '*** ERROR IN RFILRF: PINW=9999. AND IAMPLI GREATER THAN ZERO***'
            WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
            STOP
          ENDIF   !IAMPLI
          PINW=5.D0*SQRT(WAVE/(DLENG*DNPER))*(PINCEN(1)-DLENG*DNPER/2.)
          if (iampli.lt.0) pinw=pinw/(-iampli)
        ENDIF   !PINW

        IF (PINH.EQ.9999..AND.IPIN.NE.0) THEN
          IF (IAMPLI.GT.0) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
     &        '*** ERROR IN RFILRF: PINH=9999. AND IAMPLI GREATER THAN ZERO***'
            WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
            WRITE(6,*)
            WRITE(6,*)
     &        '*** ERROR IN RFILRF: PINH=9999. AND IAMPLI GREATER THAN ZERO***'
            WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
            STOP
          ENDIF   !IAMPLI
          PINH=5.D0*SQRT(WAVE/(DLENG*DNPER))*(PINCEN(1)-DLENG*DNPER/2.)
          if (iampli.lt.0) pinh=pinh/(-iampli)
        ENDIF   !PINH
      ELSE  !IUNDULATOR.NE.0
        IF (FREQLOW.EQ.-9999.) THEN
          WRITE(6,*)'*** ERROR IN RFILFR: DEFAULT FOR FREQLOW NOT ALLOWED'
          WRITE(6,*)'                     SINCE FLAG IUNDULATOR IS NOT SET'
          WRITE(6,*)'                     CHECK NAMELIST $FREQN'
          WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
          WRITE(LUNGFO,*)'*** ERROR IN RFILFR: DEFAULT FOR FREQLOW NOT ALLOWED'
          WRITE(LUNGFO,*)'                     SINCE FLAG IUNDULATOR IS NOT SET'
          WRITE(LUNGFO,*)'                     CHECK NAMELIST $FREQN'
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
          STOP
        ENDIF   !FREQLOW
        IF (FREQHIG.EQ.-9999.) THEN
          WRITE(6,*)'*** ERROR IN RFILFR: DEFAULT FOR FREQHIG NOT ALLOWED'
          WRITE(6,*)'                     SINCE FLAG IUNDULATOR IS NOT SET'
          WRITE(6,*)'                     CHECK NAMELIST $FREQN'
          WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
          WRITE(LUNGFO,*)'*** ERROR IN RFILFR: DEFAULT FOR FREQHIG NOT ALLOWED'
          WRITE(LUNGFO,*)'                     SINCE FLAG IUNDULATOR IS NOT SET'
          WRITE(LUNGFO,*)'                     CHECK NAMELIST $FREQN'
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
          STOP
        ENDIF   !FREQHIG

        IF (PINW.EQ.9999..AND.IPIN.NE.0) THEN
          WRITE(6,*)'*** ERROR IN RFILFR: DEFAULT FOR PINW NOT ALLOWED'
          WRITE(6,*)'                     SINCE FLAG IUNDULATOR IS NOT SET'
          WRITE(6,*)'                     CHECK NAMELIST $PINHOLE'
          WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
          WRITE(LUNGFO,*)'*** ERROR IN RFILFR: DEFAULT FOR PINW NOT ALLOWED'
          WRITE(LUNGFO,*)'                     SINCE FLAG IUNDULATOR IS NOT SET'
          WRITE(LUNGFO,*)'                     CHECK NAMELIST $PINHOLE'
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
          STOP
        ENDIF   !PINW

        IF (PINH.EQ.9999..AND.IPIN.NE.0) THEN
          WRITE(6,*)'*** ERROR IN RFILFR: DEFAULT FOR PINH NOT ALLOWED'
          WRITE(6,*)'                     SINCE FLAG IUNDULATOR IS NOT SET'
          WRITE(6,*)'                     CHECK NAMELIST $PINHOLE'
          WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
          WRITE(LUNGFO,*)'*** ERROR IN RFILFR: DEFAULT FOR PINH NOT ALLOWED'
          WRITE(LUNGFO,*)'                     SINCE FLAG IUNDULATOR IS NOT SET'
          WRITE(LUNGFO,*)'                     CHECK NAMELIST $PINHOLE'
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
          STOP
        ENDIF   !PINH
      ENDIF !(IUNDULATOR.NE.0)

      pinw=pinw*PINwsC
      pinh=pinh*PINHSC

C- SPECIAL CASES (FREQUENCES ARE READ FROM NAMELIST FREQN)

      IF (freqhig.le.freqlow) THEN
        IF (IFREQ2P.EQ.1) THEN
          freqhig=2.0d0*freqlow
        else IF (IFREQ2P.EQ.-1) THEN
          freqhig=1.5*freqlow
          freqlow=0.5*freqlow
         endif
      endif

      IF (IFREQ2P.EQ.0) THEN

        OPEN (UNIT=LUNFR,FILE=FILEFR,STATUS='OLD',FORM='FORMATTED',ERR=999)
c        READ(LUNFR,*,ERR=99) NFREQ

        nfreq=0
1       continue
        read(lunfr,*,end=9) freq(1)
        nfreq=nfreq+1
        goto 1
9       rewind(lunfr)

        IF(NFREQ.GT.NDFREQ) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN SR RFILFR ***'
          WRITE(LUNGFO,*)'TOO MANY FREQUENCES ON FILEFR'
          WRITE(LUNGFO,*)'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN SR RFILFR ***'
          WRITE(6,*)'TOO MANY FREQUENCES ON FILEFR'
          WRITE(6,*)'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
          STOP
        ENDIF

        DO IFR=1,NFREQ
          READ(LUNFR,*,ERR=99) DUM
          IF (IUNIT.EQ.0) THEN    ! 260194
            FREQ(IFR)=DUM
          ELSE
            FREQ(IFR)=WTOE1/DUM !260194
          ENDIF
        ENDDO

        CLOSE(LUNFR)

        nintfreq=nfreq
        freqlow=freq(1)
        freqhig=freq(nfreq)

C260194 {
        DO I=1,NFREQ
          DO J=I+1,NFREQ
            IF (IUNIT.EQ.0) THEN
              IF (FREQ(J).LT.FREQ(I)) THEN
                DUM=FREQ(J)
                FREQ(J)=FREQ(I)
                FREQ(I)=DUM
              ENDIF
            ELSE
              IF (FREQ(J).GT.FREQ(I)) THEN
                DUM=FREQ(J)
                FREQ(J)=FREQ(I)
                FREQ(I)=DUM
              ENDIF
            ENDIF
          ENDDO
        ENDDO

      else IF(IFREQ2P.GT.1) THEN

        IF(2.GT.NDFREQ) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN SR RFILFR ***'
          WRITE(LUNGFO,*)'TOO MANY FREQUENCES TO BE CALCULATED'
          WRITE(LUNGFO,*)'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN SR RFILFR ***'
          WRITE(6,*)'TOO MANY FREQUENCES TO BE CALCULATED'
          WRITE(6,*)'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
          STOP
        ENDIF

        IF (FRSCALEo.lT.0.D0) THEN
          FREQLOW=FREQLOW*frscale
          FREQHIG=FREQhig*frscale
          nintfreq=nintfreq*frscale
        else IF (FRSCALE.GT.0.D0) THEN
          FREQHIG=FREQHIG-FREQLOW
          FREQLOW=FREQLOW+FREQHIG/2.D0
          FREQLOW=FREQLOW*FRSCALE-FREQHIG/2.D0
          FREQHIG=FREQLOW+FREQHIG
        ENDIF

        NFREQ=NINTFREQ
        NFREQ0=NFREQ
        NFREQ0M=1
        NFREQ0P=NFREQ

        IF (IUNIT.EQ.0) THEN

          FREQ(1)=FREQLOW
          FREQ(2)=FREQHIG

        ELSE

          FREQ(1)=WTOE1/FREQHIG
          FREQ(2)=WTOE1/FREQLOW

        ENDIF

        IF(IFREQ2P.EQ.2)  THEN

          IF (2.D0*FREQ(1).GT.FREQ(2)) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN SR RFILFR ***'
            WRITE(LUNGFO,*)'PHOTON ENERGIES ON NAMELIST FREQN INCOMPATIBLE'
            WRITE(LUNGFO,*)'WITH FLAG IFREQ2P.'
            WRITE(LUNGFO,*)'FIRST PHOTON ENERGY MUST NOT BE HIGHER THAN'
            WRITE(LUNGFO,*)'HALF THE SECOND ONE'
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN SR RFILFR ***'
            WRITE(6,*)'PHOTON ENERGIES ON NAMELIST FREQN INCOMPATIBLE'
            WRITE(6,*)'WITH FLAG IFREQ2P.'
            WRITE(6,*)'FIRST PHOTON ENERGY MUST NOT BE HIGHER THAN'
            WRITE(6,*)'HALF THE SECOND ONE'
            STOP
          ENDIF

          NFREQ=DLOG(FREQ(2)/FREQ(1))/DLOG(2.D0)+1
          NFREQ0=NFREQ
          NFREQ0M=1
          NFREQ0P=NFREQ

        ELSE IF(IFREQ2P.GT.2) THEN

          IF (FREQ(1).GE.FREQ(2)) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN SR RFILFR ***'
            WRITE(LUNGFO,*)'PHOTON ENERGIES ON NAMELIST FREQN INCOMPATIBLE'
            WRITE(LUNGFO,*)'WITH FLAG IFREQ2P.'
            WRITE(LUNGFO,*)'FIRST PHOTON ENERGY MUST BE LOWER THAN'
            WRITE(LUNGFO,*)'THE SECOND ONE'
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN SR RFILFR ***'
            WRITE(6,*)'PHOTON ENERGIES ON NAMELIST FREQN INCOMPATIBLE'
            WRITE(6,*)'WITH FLAG IFREQ2P.'
            WRITE(6,*)'FIRST PHOTON ENERGY MUST BE LOWER THAN'
            WRITE(6,*)'THE SECOND ONE'
            STOP
          ENDIF !FREQ(1)

          NFREQ=NINTFREQ
          NFREQ0=NFREQ
          NFREQ0M=1
          NFREQ0P=NFREQ

        ENDIF   !IFREQ2P

        IF(NFREQ.GT.NDFREQ) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN SR RFILFR ***'
          WRITE(LUNGFO,*)'TOO MANY FREQUENCES TO BE CALCULATED'
          WRITE(LUNGFO,*)'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN SR RFILFR ***'
          WRITE(6,*)'TOO MANY FREQUENCES TO BE CALCULATED'
          WRITE(6,*)'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
          STOP
        ENDIF   !NFREQ

        DFREQ=(FREQ(2)-FREQ(1))/(NFREQ-1)
        DO IFR=1,NFREQ-1
          IF(IFREQ2P.EQ.2) THEN
            FREQ(IFR+1)=FREQ(IFR)*2.D0
          ELSE
            FREQ(IFR+1)=FREQ(IFR)+DFREQ
          ENDIF
        ENDDO

C260194      RETURN

      ENDIF !IFREQ2P.GT.1

      IF(iabs(IFREQ2P).EQ.1) THEN   !260194

        IF (IUNIT.EQ.0) THEN
          if (ifreq2p.eq.1) then !20220407
            FREQ(1)=FREQLOW
          else
            FREQ(1)=(FREQLOW+freqhig)/2.0d0
          endif
        ELSE
          if (ifreq2p.eq.1) then !20220407
            FREQ(1)=WTOE1/FREQLOW
          else
            FREQ(1)=(wtoe1/FREQLOW+wtoe1/freqhig)/2.0d0
          endif
        ENDIF

        nintfreq=1
        NFREQ=1
        NFREQ0=NFREQ
        NFREQ0M=1
        NFREQ0P=NFREQ

        nintfreq=nfreq
C260194  ELSE

      ENDIF !IFREQ2P

C24.3.94{
      IF (IEFOLD.NE.0.AND.IEFOLD.NE.2) THEN

        NFREQO=NFREQ
        NFREQ0=NFREQ

        IF (IFREQ2P.eq.2) THEN
          dfreq=2.0d0*ESPREAD*NSIGE*FREQlow
          IF (FREQlow-dfreq.gt.0.0d0) THEN
            freql=freqlow
            do while (freql.gt.freq(1)-dfreq)
              freql=freql/2.0d0
            enddo
            nfreq=1
            freq(1)=freql
            dfreq=2.0d0*ESPREAD*NSIGE*FREQhig
            do while (freq(nfreq).lt.freqhig+dfreq)
              nfreq=nfreq+1
              if (nfreq.gt.ndfreqp) then
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'*** ERROR IN SR RFILFR ***'
                WRITE(LUNGFO,*)'TOO MANY FREQUENCES TO BE CALCULATED'
                WRITE(LUNGFO,*)'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
                WRITE(6,*)
                WRITE(6,*)'*** ERROR IN SR RFILFR ***'
                WRITE(6,*)'TOO MANY FREQUENCES TO BE CALCULATED'
                WRITE(6,*)'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
                STOP
              endif
              freq(nfreq)=freq(nfreq-1)*2.0d0
            enddo
          else
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN RFILRF ***'
            WRITE(LUNGFO,*)
     &        'NEGATIVE OR ZERO PHOTON ENERGY OCCURED WHILE EXTENDING ENERGY RANGE DUE TO FLAG IEFOLD'
            WRITE(LUNGFO,*)'CHECK INPUT FILE, INCREASE FREQLOW OR NINTFREQ'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN RFILRF ***'
            WRITE(6,*)
     &        'NEGATIVE OR ZERO PHOTON ENERGY OCCURED WHILE EXTENDING ENERGY RANGE DUE TO FLAG IEFOLD'
            WRITE(6,*)'CHECK INPUT FILE, INCREASE FREQLOW OR NINTFREQ'
            WRITE(6,*)
            STOP
          ENDIF
        else IF (IFREQ2P.eq.1.or.freqlow.eq.freqhig) THEN
          dfreq=2.0d0*ESPREAD*NSIGE*FREQlow
          freq(1)=freqlow-dfreq
          dfreq=dfreq/5.0d0
          nfreq=11
          do ifreq=2,nfreq
            freq(ifreq)=freq(ifreq-1)+dfreq
          enddo
        else IF (IFREQ2P.eq.-1) THEN
          dfreq=2.0d0*ESPREAD*NSIGE*(FREQlow+freqhig)/2.
          freq(1)=(freqlow+freqhig)/2.-dfreq
          dfreq=dfreq/5.0d0
          nfreq=11
          do ifreq=2,nfreq
            freq(ifreq)=freq(ifreq-1)+dfreq
          enddo
        else !ifreq2p

          DFREQ=FREQ(nfreq)-FREQ(nfreq-1)
          DEXTEND=2.0d0*ESPREAD*NSIGE*FREQ(NFREQ)
          INCREAS=DMAX1(2.D0,DEXTEND/DFREQ+1.D0)
          NFREQEP=INCREAS
          NFREQ=NFREQ+INCREAS

          IF(NFREQ.GT.NDFREQ) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN SR RFILFR ***'
            WRITE(LUNGFO,*)'TOO MANY FREQUENCES TO BE CALCULATED'
            WRITE(LUNGFO,*)'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN SR RFILFR ***'
            WRITE(6,*)'TOO MANY FREQUENCES TO BE CALCULATED'
            WRITE(6,*)'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
            STOP
          ENDIF   !NFREQ

          DO IFREQ=1,INCREAS
            FREQ(NFREQO+IFREQ)=FREQ(NFREQO)+DFREQ*IFREQ
          ENDDO

          DFREQ=FREQ(2)-FREQ(1)
          DEXTEND=2.0d0*ESPREAD*NSIGE*FREQ(1)
          INCREAS=DMAX1(2.D0,DEXTEND/DFREQ+1.D0)
          NFREQEM=INCREAS
          NFREQ0M=NFREQEM+1
          NFREQ0P=NFREQ0M+NFREQ0-1
          NFREQO=NFREQ
          NFREQ=NFREQ+INCREAS

          IF(NFREQ.GT.NDFREQ) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN SR RFILFR ***'
            WRITE(LUNGFO,*)'TOO MANY FREQUENCES TO BE CALCULATED'
            WRITE(LUNGFO,*)'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN SR RFILFR ***'
            WRITE(6,*)'TOO MANY FREQUENCES TO BE CALCULATED'
            WRITE(6,*)'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
            STOP
          ENDIF   !NFREQ

          DO IFREQ=NFREQO,1,-1
            FREQ(IFREQ+INCREAS)=FREQ(IFREQ)
          ENDDO
          DO IFREQ=1,INCREAS
            FREQ(INCREAS+1-IFREQ)=FREQ(INCREAS+1)-DFREQ*IFREQ
            IF (FREQ(INCREAS+1-IFREQ).LE.0.0) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** ERROR IN RFILRF ***'
              WRITE(LUNGFO,*)
     &          'NEGATIVE OR ZERO PHOTON ENERGY OCCURED WHILE EXTENDING ENERGY RANGE DUE TO FLAG IEFOLD'
              WRITE(LUNGFO,*)'CHECK INPUT FILE'
              WRITE(LUNGFO,*)
              WRITE(6,*)
              WRITE(6,*)'*** ERROR IN RFILRF ***'
              WRITE(6,*)
     &          'NEGATIVE OR ZERO PHOTON ENERGY OCCURED WHILE EXTENDING ENERGY RANGE DUE TO FLAG IEFOLD'
              WRITE(6,*)'CHECK INPUT FILE'
              WRITE(6,*)
              STOP
            ENDIF
          ENDDO
        ENDIF !ifreq2p

      ENDIF   !IEFOLD
C24.3.94}

      DO IFR=1,NFREQ
        WELLEN(IFR)=WTOE1/FREQ(IFR)
      ENDDO

C260194 }


      RETURN

999   CONTINUE

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'*** ERROR IN SR RFILFR ***'
      WRITE(LUNGFO,*)'FILE OPENING ERROR '
      WRITE(LUNGFO,*)'FILE, UNIT:'
      WRITE(LUNGFO,*)FILEFR, LUNFR
      WRITE(6,*)
      WRITE(6,*)'*** ERROR IN SR RFILFR ***'
      WRITE(6,*)'FILE OPENING ERROR '
      WRITE(6,*)'FILE, UNIT:'
      WRITE(6,*)FILEFR, LUNFR
      STOP
99    CONTINUE
      WRITE(6,*)
      WRITE(6,*)'*** ERROR IN SR RFILFR ***'
      WRITE(6,*)'FILE READING ERROR '
      WRITE(6,*)'FILE, UNIT:'
      WRITE(6,*)FILEFR, LUNFR
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'*** ERROR IN SR RFILFR ***'
      WRITE(LUNGFO,*)'FILE READING ERROR '
      WRITE(LUNGFO,*)'FILE, UNIT:'
      WRITE(LUNGFO,*)FILEFR, LUNFR
      STOP
      END
