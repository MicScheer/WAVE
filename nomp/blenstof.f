*CMZ :  3.07/00 15/03/2019  13.20.05  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.35/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.34/09 26/09/2001  12.12.32  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  12.20.49  by  Michael Scheer
*CMZ :  2.16/08 24/10/2000  14.28.41  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.08.20  by  Michael Scheer
*CMZ :  2.13/03 12/01/2000  16.31.33  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.04.42  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.51.13  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.58  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.05  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BLENSTOF(ISTOK,IFREQ)
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

      DOUBLE PRECISION DSUM,RPHI,SUMP

      INTEGER IFREQ,ISTOK,IY,IZ,IOBSV,IIY,IIZ
     &  ,ICAL,IR,IP,MR,MP,KDUM,IERR

      DOUBLE PRECISION SUMZ(NDOBSVYP),S2(NDOBSVYP),SUM,SUMY(NDOBSVZP)
      DOUBLE PRECISION OBSVYF(NDOBSVYP),DIA

      DOUBLE PRECISION R(NDOBSVZP),PHI(NDOBSVYP)
     &  ,FPHI(NDOBSVYP)
     &  ,SZ(NDOBSVZP),SY(NDOBSVYP)
     &  ,FZ(NDOBSVZP),FY(NDOBSVYP),DPHI,DR,X,Y
      DOUBLE PRECISION OBSVZF(NDOBSVZP)

      DATA ICAL/0/

C--- TAKE INNER EDGE OF PINHOLE INTO ACCOUNT, I.E. SET MOBSVZ,MOBVY,MOBSV
C    TO ORIGINAL VALUES. THEY HAVE BEEN OVERWRITTEN IN SR WFOLINT

      MOBSVZ=MOBSVZ-2*MMEDGEZ
      MOBSVY=MOBSVY-2*MMEDGEY
      MOBSV=MOBSVZ*MOBSVY

      IF (IPINCIRC.EQ.0) THEN

        IF (IF1DIM.NE.2) THEN

C--- INTEGRATION ALONG HORIZONTAL AXIS Z

          IIY=0
          DO IY=(NOBSVY-MOBSVY-2*MMEDGEY)/2+1,
     &        (NOBSVY-MOBSVY-2*MMEDGEY)/2+MOBSVY+2*MMEDGEY

            IIY=IIY+1
            OBSVYF(IIY)=OBSVY(IY)

            IF(MOBSVZ.GT.1) THEN

              SUMZ(IIY)=0.0d0
              DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ-1

                IOBSV=(IY-1)*NOBSVZ+IZ
                IOBFR=IOBSV+NOBSV*(IFREQ-1)
                DSUM=
     &            OBSVDZ*0.5D0
     &            *(STOKESF(ISTOK,IOBFR)+STOKESF(ISTOK,IOBFR+1))
     &            -OBSVDZ**3/24.D0
     &            *(SPCOEFM(IOBSV)
     &            + SPCOEFM(IOBSV+1))

                SUMZ(IIY)=SUMZ(IIY)+DSUM

              ENDDO   !IZ

            ELSE

              IOBSV=(IY-1)*NOBSVZ+(NOBSVZ-MOBSVZ)/2+1
              SUMZ(IIY)=OBSVDZ*STOKESF(ISTOK,IOBSV+NOBSV*(IFREQ-1))
            ENDIF

          ENDDO !IY

        ELSE   !IF1DIM.EQ.2

C--- INTEGRATION ALONG HORIZONTAL AXIS Z

          IIY=0
          DO IY=(NOBSVY-MOBSVY-2*MMEDGEY)/2+1,
     &        (NOBSVY-MOBSVY-2*MMEDGEY)/2+MOBSVY+2*MMEDGEY


            IIY=IIY+1
            OBSVYF(IIY)=OBSVY(IY)

            IOBSV=(IY-1)*NOBSVZ+(NOBSVZ-MOBSVZ)/2+1
            DIA=ABS((PINR-(PINCEN(2)-OBSV(2,IOBSV)))
     &        *(PINR+(PINCEN(2)-OBSV(2,IOBSV))))
            IF (DIA.GT.0.D0) THEN
              DIA=2.D0*SQRT(DIA)
            ELSE
              DIA=0.D0
            ENDIF
            SUMZ(IIY)=DIA*STOKESF(ISTOK,IOBSV+NOBSV*(IFREQ-1))

          ENDDO !IY

        ENDIF   !IF1DIM.EQ.2

C--- INTEGRATION ALONG VERTICAL AXIS Y

        CALL FSPLINDX(OBSVDY,SUMZ,MOBSVY+2*MMEDGEY,0.D0,0.D0,S2)

        IF(MOBSVY.GT.1) THEN

          SUM=0.0d0

          DO IY=MMEDGEY+1,MMEDGEY+MOBSVY-1

            DSUM=
     &        OBSVDY*0.5D0
     &        *(SUMZ(IY)+SUMZ(IY+1))
     &        -OBSVDY**3/24.D0
     &        *(S2(IY)+S2(IY+1))

            SUM=SUM+DSUM

          ENDDO

        ELSE IF (IF1DIM.EQ.2) THEN

          SUM=PI1*PINR*SUMZ(MMEDGEY+1)/2.D0

        ELSE

          SUM=OBSVDY*SUMZ(MMEDGEY+1)

        ENDIF

      ELSE  !IPINCIRC

C--- INTEGRATION OVER PHI

        IF (IRPHI.NE.0) THEN !INTEGRATION WITH RESPECT TO POLAR COORDINATES

          IF (ICAL.EQ.0) THEN

            DR=DMIN1(OBSVDZ,OBSVDY)
            MR=NINT(PINR/DR)+1
            DR=PINR/(MR-1)
            MEDGER=MIN(MMEDGEZ,MMEDGEY)
            MP=MOBSVY

            IF (MR.LT.1) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** ERROR IN BLENSTOF ***'
              WRITE(LUNGFO,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
              WRITE(LUNGFO,*)'INCREASE PARAMETER MPINZ IN NAMELIST PINHOLE'
              WRITE(LUNGFO,*)
              WRITE(6,*)
              WRITE(6,*)'*** ERROR IN BLENSTOF ***'
              WRITE(6,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
              WRITE(6,*)'INCREASE PARAMETER MPINZ IN NAMELIST PINHOLE'
              WRITE(6,*)
            ENDIF
            IF (MP.LT.4) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** ERROR IN BLENSTOF ***'
              WRITE(LUNGFO,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
              WRITE(LUNGFO,*)'INCREASE PARAMETER MPINY IN NAMELIST PINHOLE'
              WRITE(LUNGFO,*)
              WRITE(6,*)
              WRITE(6,*)'*** ERROR IN BLENSTOF ***'
              WRITE(6,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
              WRITE(6,*)'INCREASE PARAMETER MPINY IN NAMELIST PINHOLE'
              WRITE(6,*)
            ENDIF

            DPHI=2.D0*PI1/(MP-1)
            DO IP=1,MP
              PHI(IP)=(IP-1)*DPHI
            ENDDO

            DO IR=1,MR+MEDGER
              R(IR)=(IR-1)*DR
            ENDDO

            DO IR=2,MR+MEDGER
              DO IP=1,MP
                XC(IP+(IR-1)*NOBSVY)=R(IR)*DCOS(PHI(IP))+PINCEN(3)
                YC(IP+(IR-1)*NOBSVY)=R(IR)*DSIN(PHI(IP))+PINCEN(2)
              ENDDO !IP
            ENDDO !IR

            ICAL=1

          ENDIF !ICAL

C--- INTERPOLATION OF INTENSITY FOR CIRCULAR GRID

          DO IR=2,MR+MEDGER
            DO IP=1,MP

              X=XC(IP+(IR-1)*NOBSVY)
              Y=YC(IP+(IR-1)*NOBSVY)

              IIY=0
              DO IY=(NOBSVY-MOBSVY-2*MMEDGEY)/2+1,
     &            (NOBSVY-MOBSVY-2*MMEDGEY)/2+MOBSVY+2*MMEDGEY
                IIY=IIY+1

                IIZ=0
                DO IZ=(NOBSVZ-MOBSVZ-2*MMEDGEZ)/2+1,
     &              (NOBSVZ-MOBSVZ-2*MMEDGEZ)/2+MOBSVZ+2*MMEDGEZ
                  IIZ=IIZ+1
                  IOBSV=(IY-1)*NOBSVZ+IZ
                  FZ(IIZ)=STOKESF(ISTOK,IOBSV+NOBSV*(IFREQ-1))
                  SZ(IIZ)=SPCOEFM(IOBSV)
                  OBSVZF(IIZ)=OBSVZ(IZ)
                ENDDO !IZ

                CALL SPLINZY(IIZ,X,FY(IIY),OBSVZF,FZ,SZ,KDUM)

                OBSVYF(IIY)=OBSVY(IY)

              ENDDO !IY

              CALL FSPLINDX(OBSVDY,FY,IIY,0.D0,0.D0,SY)
              CALL SPLINZY(IIY,Y,FPHIR(IP+(IR-1)*NOBSVY),OBSVYF,FY,SY,KDUM)

            ENDDO !IP
          ENDDO !IR

C--- DO THE INTEGRATION OF FPHIR OVER PHI AND R

          SUM=0.0D0
          SUMY(1)=0.0D0
          DO IR=2,MR+MEDGER !FIRST RADIUS IS ZERO

            DO IP=1,MP
              FPHI(IP)=FPHIR(IP+(IR-1)*NOBSVY)
            ENDDO   !IP

            CALL FSPLPER(DPHI,FPHI,MP,SY)

            SUMY(IR)=0.0D0
            RPHI=R(IR)*DPHI
            DO IP=1,MP-1

              DSUM=
     &          RPHI*0.5D0*(FPHI(IP)+FPHI(IP+1))
     &          -RPHI**3/24.D0*(SY(IP)+SY(IP+1))

              SUMY(IR)=SUMY(IR)+DSUM

            ENDDO   !IP

          ENDDO !IR

          CALL FSPLINDX(DR,SUMY,MR+MEDGER,0.D0,0.D0,SZ)

          SUM=0.0d0
          DO IR=1,MR-1

            DSUM=
     &        DR*0.5D0
     &        *(SUMY(IR)+SUMY(IR+1))
     &        -DR**3/24.D0
     &        *(SZ(IR)+SZ(IR+1))

            SUM=SUM+DSUM

          ENDDO

        ELSE  !IRPHI

          DO IOBSV=1,NOBSV
            FPHIR(IOBSV)=STOKESF(ISTOK,IOBSV+NOBSV*(IFREQ-1))
          ENDDO !IOBSV

          CALL CIRCPIN(NOBSVZ,NOBSVY,MOBSVZ,MOBSVY,FPHIR,SUM,SUMP,-ISTOK,IFREQ,IERR)

        ENDIF !IRPHI

      ENDIF !PINCIRC

      WSTOKESF(ISTOK,IFREQ)=SUM

      RETURN
      END
