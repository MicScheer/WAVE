*CMZ :  3.08/01 03/04/2019  14.39.38  by  Michael Scheer
*CMZ :  3.07/00 15/03/2019  15.22.31  by  Michael Scheer
*CMZ :  3.03/02 27/02/2017  13.51.45  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  16.22.59  by  Michael Scheer
*CMZ :  2.52/09 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.35/01 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.34/09 24/09/2001  12.08.47  by  Michael Scheer
*CMZ :  2.34/01 25/06/2001  14.39.42  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  12.14.41  by  Michael Scheer
*CMZ :  2.31/00 23/04/2001  18.27.11  by  Michael Scheer
*CMZ :  2.20/09 23/03/2001  11.01.07  by  Michael Scheer
*CMZ :  2.16/08 24/10/2000  14.08.07  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.08.20  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.22.27  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.04.41  by  Michael Scheer
*CMZ :  1.00/00 31/07/97  17.36.12  by  Michael Scheer
*CMZ : 00.02/05 18/03/97  15.48.44  by  Michael Scheer
*CMZ : 00.02/04 26/02/97  12.07.43  by  Michael Scheer
*CMZ : 00.02/00 10/12/96  18.07.52  by  Michael Scheer
*CMZ : 00.01/09 01/09/95  12.58.16  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.49.44  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.38  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.04  by  Michael Scheer
*-- Author :
      SUBROUTINE BLENDE(ISOUR,kfreq)

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

C--- INTEGRATES THE SPLINES THAT INTERPOLATE THE INTENSITY INSIDE THE PINHOLE

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

      INTEGER kfreq,ISOUR,IY,IZ,IOBSV,IIY,ICAL,IR,IP,NR,MR,MP,KDUM,IERR
      INTEGER JCAL,IWBLEN,IDUM
      INTEGER IWRPHIS,IWRPHIF,IWSOUR,IWFREQ,ILIOBFR1

      DOUBLE PRECISION DSUM,RPHI,DIA
      DOUBLE PRECISION SUMZ(NDOBSVYP),S2(NDOBSVYP),SUM,OBSVYF(NDOBSVYP)
      DOUBLE PRECISION SUMY(NDOBSVZP)
      DOUBLE PRECISION SUMZP(NDOBSVYP),S2P(NDOBSVYP),SUMP,DSUMP
      DOUBLE PRECISION R(NDOBSVZP),PHI(NDOBSVYP)
     &  ,FPHI(NDOBSVYP)
     &  ,SZ(NDOBSVZP),SY(NDOBSVYP)
     &  ,FZ(NDOBSVZP),FY(NDOBSVYP),DPHI,DR,X,Y

      DATA ICAL/0/,JCAL/0/,iwblen/0/

      if (ipin.eq.3) then
        ILIOBFR=ISOUR+NSOURCE*(kfreq-1)
        if (ipincirc.eq.0) then
          WFLUX(ILIOBFR)=SPEC(ILIOBFR)*pinw*pinh
        else
          WFLUX(ILIOBFR)=SPEC(ILIOBFR)*pinr**2*pi1
        endif
        return
      endif

c      IWBLEN=0

      IF (IPINCIRC.EQ.0) THEN

        IF (IF1DIM.NE.2) THEN

C--- INTEGRATION ALONG HORIZONTAL AXIS Z

          IIY=0
C290693  DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY
          DO IY=1,NOBSVY

            IIY=IIY+1
            SUMZ(IIY)=0.0
            SUMZP(IIY)=0.0
            OBSVYF(IIY)=OBSVY(IY)

            IF(MOBSVZ.GT.1) THEN

              DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ-1

                IOBSV=(IY-1)*NOBSVZ+IZ

                ILIOBFR=ISOUR+NSOURCE*(IOBSV-1+NOBSV*(kfreq-1))
                ILIOBFR1=ISOUR+NSOURCE*(IOBSV+NOBSV*(kfreq-1))

                DSUM=OBSVDZ*0.5D0
     &            *(SPEC(ILIOBFR)+SPEC(ILIOBFR1))
     &            -OBSVDZ**3/24.D0*(SPCOEF(IOBSV)+SPCOEF(IOBSV+1))

                IF(
     &              (IWSOUR.NE.ISOUR.OR.IWFREQ.NE.kfreq)
     &              .AND.
     &              DSUM.LT.0.0) THEN

                  IF (IWBLEN.EQ.0) THEN
                    WRITE(LUNGFO,*)
                    WRITE(LUNGFO,*)
                    WRITE(LUNGFO,*)'*** WARNING IN BLENDE ***'
                    WRITE(LUNGFO,*)
     &                'SPLINE INTEGRATION FAILED, RESULTS NOT RELIABLE'
                    WRITE(LUNGFO,*)
                    WRITE(LUNGFO,*)
                  ENDIF
                  IWSOUR=ISOUR
                  IWFREQ=kfreq
                  IW_BLEN=1
                  IWBLEN=1
                  DO IDUM=1,IIY
                    SUMZP(IDUM)=SUMZ(IDUM)
                  ENDDO
                ENDIF !IWSOUR

                IF (IWBLEN.NE.0) SUMZP(IIY)=SUMZP(IIY)+DABS(DSUM)
                SUMZ(IIY)=SUMZ(IIY)+DSUM

              ENDDO   !IZ

            ELSE  !MOBSVZ

              IOBSV=(IY-1)*NOBSVZ+(NOBSVZ-MOBSVZ)/2+1
              SUMZ(IIY)=OBSVDZ*SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(kfreq-1)))

              IF (IWBLEN.NE.0) SUMZP(IIY)=SUMZ(IIY)

            ENDIF !MOBSVZ

          ENDDO !IY

C--- INTEGRATION ALONG VERTICAL AXIS Y

          CALL FSPLINDX(OBSVDY,SUMZ,NOBSVY,0.D0,0.D0,S2)
          IF (IWBLEN.NE.0)
     &      CALL FSPLINDX(OBSVDY,SUMZP,NOBSVY,0.D0,0.D0,S2P)

          IF(MOBSVY.GT.1) THEN

            SUM=0.0
            SUMP=0.0
            DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY-1

              DSUM=
     &          OBSVDY*0.5D0
     &          *(SUMZ(IY)+SUMZ(IY+1))
     &          -OBSVDY**3/24.D0
     &          *(S2(IY)+S2(IY+1))

              IF (IWBLEN.NE.0) THEN
                DSUMP=
     &            OBSVDY*0.5D0
     &            *(SUMZP(IY)+SUMZP(IY+1))
     &            -OBSVDY**3/24.D0
     &            *(S2P(IY)+S2P(IY+1))
              ENDIF

              IF(
     &            (IWSOUR.NE.ISOUR.OR.IWFREQ.NE.kfreq)
     &            .AND.
     &            DSUM.LT.0.0) THEN
                IF (IWBLEN.EQ.0) THEN
                  WRITE(LUNGFO,*)
                  WRITE(LUNGFO,*)
                  WRITE(LUNGFO,*)'*** WARNING IN BLENDE ***'
                  WRITE(LUNGFO,*)
     &              'SPLINE INTEGRATION FAILED, RESULTS NOT RELIABLE'
                  WRITE(LUNGFO,*)
                  WRITE(LUNGFO,*)
                ENDIF

                IWSOUR=ISOUR
                IWFREQ=kfreq
                IW_BLEN=1
                IWBLEN=1

                DO IDUM=1,NOBSVY
                  S2P(IDUM)=S2(IDUM)
                  SUMZP(IDUM)=SUMZ(IDUM)
                  SUMP=SUM
                  DSUMP=DSUM
                ENDDO

              ENDIF !IWSOUR

              IF (IWBLEN.NE.0) SUMP=SUMP+DABS(DSUMP)

              SUM=SUM+DSUM

            ENDDO !IY

          ELSE !MOBSVY

            SUM=OBSVDY*SUMZ(NOBSVY/2+1)

            IF (IWBLEN.NE.0) SUMP=OBSVDY*SUMZP(NOBSVY/2+1)

          ENDIF !MOBSVY

        ELSE !(IF1DIM.EQ.2)

C--- INTEGRATION ALONG HORIZONTAL AXIS Z

          DO IY=1,NOBSVY

            SUMZP(IY)=0.0
            OBSVYF(IY)=OBSVY(IY)

            DIA=ABS((PINR-(PINCEN(2)-OBSV(2,IY)))
     &        *(PINR+(PINCEN(2)-OBSV(2,IY))))
            IF (DIA.GT.0.D0) THEN
              DIA=2.D0*SQRT(DIA)
            ELSE
              DIA=0.D0
            ENDIF
            SUMZ(IY)=DIA*SPEC(ISOUR+NSOURCE*(IY-1+NOBSV*(kfreq-1)))

            IF (IWBLEN.NE.0) SUMZP(IY)=SUMZ(IY)

          ENDDO !IY

        ENDIF  !IF1DIM

C--- INTEGRATION ALONG VERTICAL AXIS Y

        CALL FSPLINDX(OBSVDY,SUMZ,NOBSVY,0.D0,0.D0,S2)
        IF (IWBLEN.NE.0)
     &    CALL FSPLINDX(OBSVDY,SUMZP,NOBSVY,0.D0,0.D0,S2P)

        IF(MOBSVY.GT.1) THEN

          SUM=0.0
          SUMP=0.0
          DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY-1


            DSUM=
     &        OBSVDY*0.5D0
     &        *(SUMZ(IY)+SUMZ(IY+1))
     &        -OBSVDY**3/24.D0
     &        *(S2(IY)+S2(IY+1))

            IF (IWBLEN.NE.0) THEN
              DSUMP=
     &          OBSVDY*0.5D0
     &          *(SUMZP(IY)+SUMZP(IY+1))
     &          -OBSVDY**3/24.D0
     &          *(S2P(IY)+S2P(IY+1))
            ENDIF

            IF(
     &          (IWSOUR.NE.ISOUR.OR.IWFREQ.NE.kfreq)
     &          .AND.
     &          DSUM.LT.0.0) THEN
              IF (IWBLEN.EQ.0) THEN
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'*** WARNING IN BLENDE ***'
                WRITE(LUNGFO,*)
     &            'SPLINE INTEGRATION FAILED, RESULTS NOT RELIABLE'
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)
              ENDIF

              IWSOUR=ISOUR
              IWFREQ=kfreq
              IW_BLEN=1
              IWBLEN=1
              DO IDUM=1,NOBSVY
                S2P(IDUM)=S2(IDUM)
                SUMZP(IDUM)=SUMZ(IDUM)
                SUMP=SUM
                DSUMP=DSUM
              ENDDO

            ENDIF !IWSOUR

            IF (IWBLEN.NE.0) SUMP=SUMP+DABS(DSUMP)

            SUM=SUM+DSUM

          ENDDO !IY

        ELSE IF (IF1DIM.EQ.2) THEN

          SUM=PI1*PINR*SUMZ(NOBSVY/2+1)/2.D0

          IF (IWBLEN.NE.0) SUMP=PI1*PINR*SUMZP(NOBSVY/2+1)/2.D0

        ELSE !MOBSVY

          SUM=OBSVDY*SUMZ(NOBSVY/2+1)

          IF (IWBLEN.NE.0) SUMP=OBSVDY*SUMZP(NOBSVY/2+1)

        ENDIF !MOBSVY

      ELSE  !IPINCIRC (BZW. IF1DIM.EQ.1)

        IF (IRPHI.NE.0) THEN !INTEGRATION WITH RESPECT TO POLAR COORDINATES

C--- INTEGRATION OVER PHI

          IF (ICAL.EQ.0) THEN

            DR=DMIN1(OBSVDZ,OBSVDY)
            MR=NINT(PINR/DR)+1
            DR=PINR/(MR-1)
            MEDGER=MIN( MEDGEZ, MEDGEY)
            NR=MR+MEDGER
            MP=MOBSVY

            IF (MR.LT.1) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** ERROR IN BLENDE ***'
              WRITE(LUNGFO,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
              WRITE(LUNGFO,*)'INCREASE PARAMETER MPINZ IN NAMELIST PINHOLE'
              WRITE(LUNGFO,*)
              WRITE(6,*)
              WRITE(6,*)'*** ERROR IN BLENDE ***'
              WRITE(6,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
              WRITE(6,*)'INCREASE PARAMETER MPINZ IN NAMELIST PINHOLE'
              WRITE(6,*)
            ENDIF
            IF (MP.LT.4) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** ERROR IN BLENDE ***'
              WRITE(LUNGFO,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
              WRITE(LUNGFO,*)'INCREASE PARAMETER MPINY IN NAMELIST PINHOLE'
              WRITE(LUNGFO,*)
              WRITE(6,*)
              WRITE(6,*)'*** ERROR IN BLENDE ***'
              WRITE(6,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
              WRITE(6,*)'INCREASE PARAMETER MPINY IN NAMELIST PINHOLE'
              WRITE(6,*)
            ENDIF

            DPHI=2.D0*PI1/(MP-1)
            DO IP=1,MP
              PHI(IP)=(IP-1)*DPHI
            ENDDO

            DO IR=1,NR
              R(IR)=(IR-1)*DR
            ENDDO

            DO IR=2,NR
              DO IP=1,MP
                XC(IP+(IR-1)*NOBSVY)=R(IR)*DCOS(PHI(IP))+PINCEN(3)
                YC(IP+(IR-1)*NOBSVY)=R(IR)*DSIN(PHI(IP))+PINCEN(2)
              ENDDO !IP
            ENDDO !IR

            ICAL=1

          ENDIF !ICAL

C--- INTERPOLATION OF INTENSITY FOR CIRCULAR GRID

          DO IR=2,NR
            DO IP=1,MP

              X=XC(IP+(IR-1)*NOBSVY)
              Y=YC(IP+(IR-1)*NOBSVY)

              DO IY=1,NOBSVY

                DO IZ=1,NOBSVZ
                  IOBSV=(IY-1)*NOBSVZ+IZ
                  FZ(IZ)=SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(kfreq-1)))
                  SZ(IZ)=SPCOEF(IOBSV)
                ENDDO !IZ

                CALL SPLINZY(NOBSVZ,X,FY(IY),OBSVZ,FZ,SZ,KDUM)

              ENDDO !IY

              CALL FSPLINDX(OBSVDY,FY,NOBSVY,0.D0,0.D0,SY)
              CALL SPLINZY(NOBSVY,Y,FPHIR(IP+(IR-1)*NOBSVY),OBSVY,FY,SY,KDUM)

              IF(
     &          (IWRPHIS.NE.ISOUR.OR.IWRPHIF.NE.kfreq)
     &          .AND.
     &          FPHIR(IP+(IR-1)*NOBSVY).LT.0.0) THEN
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'*** WARNING IN BLENDE ***'
                WRITE(LUNGFO,*)
     &            'SPLINE INTERPOLATION FOR OPTION IRPHI FAILED, RESULTS NOT RELIABLE'
                WRITE(LUNGFO,*)'SOURCE POINT AND PHOTON ENERGY:'
     &            ,ISOUR,SNGL(FREQ(kfreq))
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)
                IWRPHIS=ISOUR
                IWRPHIF=kfreq
                IW_BLEN=1
                IWBLEN=1
              ENDIF !IWPHIR

            ENDDO !IP
          ENDDO !IR

C--- DO THE INTEGRATION OF FPHIR OVER PHI AND R

          SUM=0.D0
          SUMY(1)=0.0
          DO IR=2,NR  !FIRST RADIUS IS ZERO

            DO IP=1,MP
              FPHI(IP)=FPHIR(IP+(IR-1)*NOBSVY)
            ENDDO   !IP

            CALL FSPLPER(DPHI,FPHI,MP,SY)

            SUMY(IR)=0.D0
            RPHI=R(IR)*DPHI
            DO IP=1,MP-1

              DSUM=
     &          RPHI*0.5D0*(FPHI(IP)+FPHI(IP+1))
     &          -RPHI**3/24.D0*(SY(IP)+SY(IP+1))

              IF(
     &            (IWSOUR.NE.ISOUR.OR.IWFREQ.NE.kfreq)
     &            .AND.
     &            DSUM.LT.0.0) THEN
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'*** WARNING IN BLENDE ***'
                WRITE(LUNGFO,*)
     &            'SPLINE INTEGRATION FAILED, RESULTS NOT RELIABLE'
                WRITE(LUNGFO,*)'SOURCE POINT AND PHOTON ENERGY:'
     &            ,ISOUR,SNGL(FREQ(kfreq))
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)
                IWSOUR=ISOUR
                IWFREQ=kfreq
                IW_BLEN=1
                IWBLEN=1
              ENDIF !IWSOUR

              SUMY(IR)=SUMY(IR)+DSUM

            ENDDO   !IP

          ENDDO !IR

          CALL FSPLINDX(DR,SUMY,NR,0.D0,0.D0,SZ)

          SUM=0.0
          DO IR=1,MR-1

            DSUM=
     &        DR*0.5D0
     &        *(SUMY(IR)+SUMY(IR+1))
     &        -DR**3/24.D0
     &        *(SZ(IR)+SZ(IR+1))

            IF(
     &          (IWSOUR.NE.ISOUR.OR.IWFREQ.NE.kfreq)
     &          .AND.
     &          DSUM.LT.0.0) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** WARNING IN BLENDE ***'
              WRITE(LUNGFO,*)
     &          'SPLINE INTEGRATION FAILED, RESULTS NOT RELIABLE'
              WRITE(LUNGFO,*)'SOURCE POINT AND PHOTON ENERGY:'
     &          ,ISOUR,SNGL(FREQ(kfreq))
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              IWSOUR=ISOUR
              IWFREQ=kfreq
              IW_BLEN=1
              IWBLEN=1
            ENDIF !IWSOUR

            SUM=SUM+DSUM

          ENDDO

        ELSE  !IRPHI

          DO IOBSV=1,NOBSV
            FPHIR(IOBSV)=SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(kfreq-1)))
          ENDDO !IOBSV

          CALL CIRCPIN(NOBSVZ,NOBSVY,MOBSVZ,MOBSVY,FPHIR,SUM,SUMP,ISOUR,kfreq,
     &      IERR)

          IF (IERR.NE.0) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'SR CIRCPIN HAS BEEN CALLED BY SR BLENDE WITH ERRORS'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'SR CIRCPIN HAD BEEN CALLED BY SR BLENDE WITH ERRORS'
            WRITE(6,*)
          ENDIF

        ENDIF !IRPHI

      ENDIF !IPINCIRC

      WFLUX(ISOUR+NSOURCE*(kfreq-1))=SUM

      IF (IRPHI.EQ.0.AND.(IWBLEN.NE.0.OR.IERR.NE.0)) THEN
        IF (JCAL.EQ.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'      *** SUBROUTINE BLENDE:'
          WRITE(LUNGFO,*)
     &      '      LINES INDICATED BY * SHOW A RAW ESTIMATE OF ERRORS DUE TO'
          WRITE(LUNGFO,*)
     &      '      SPLINE FAILURE IF REL. ERROR .GT. 1E-5 (FIRST NUMBER IS SOURCE)'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'      *** SUBROUTINE BLENDE:'
          WRITE(6,*)
     &      '      LINES INDICATED BY * SHOW A RAW ESTIMATE OF ERRORS DUE TO'
          WRITE(6,*)
     &      '      SPLINE FAILURE IF REL. ERROR .GT. 1E-5 (FIRST NUMBER IS SOURCE)'
          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*)
     &      '      source, energy, flux, flux+error, ratio:'
          JCAL=1
        ENDIF   !JCAL
        IF (SUMP.NE.0.D0) THEN
          DSUM=SUM/SUMP
        ELSE
          DSUM=-9999.
        ENDIF
        IF (DABS(DSUM-1.D0).GT.1.D-5) THEN
          WRITE(LUNGFO,*)'*',ISOUR,
     &      SNGL(FREQ(kfreq)),SNGL(SUM),SNGL(SUMP),SNGL(DSUM)
          WRITE(6,*)'*',ISOUR,
     &      SNGL(FREQ(kfreq)),SNGL(SUM),SNGL(SUMP),SNGL(DSUM)
        ENDIF
      ENDIF !IWBLEN

      RETURN
      END
