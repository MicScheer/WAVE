*CMZ :  3.08/01 04/04/2019  12.00.57  by  Michael Scheer
*CMZ :  3.07/00 16/03/2019  15.27.40  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.35/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.34/09 24/09/2001  12.09.00  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  12.20.12  by  Michael Scheer
*CMZ :  2.16/08 24/10/2000  14.08.42  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  16.26.19  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.08.20  by  Michael Scheer
*CMZ :  2.13/04 21/01/2000  12.17.14  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.22.27  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.04.41  by  Michael Scheer
*CMZ : 00.02/04 25/02/97  17.36.07  by  Michael Scheer
*CMZ : 00.02/00 10/12/96  18.09.03  by  Michael Scheer
*CMZ : 00.01/09 01/09/95  13.01.10  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.50.50  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.44  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.06  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE blende_omp(ISOUR,kfreq)

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

      use circpinmod

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


      INTEGER IWRPHIS,IWRPHIF,IWSOUR,IWFREQ
      DOUBLE PRECISION DSUM,RPHI

      INTEGER kfreq,ISOUR,IY,IZ,IOBSV,IIY,IIZ
     &  ,ICAL,IR,IP,MR,MP,KDUM,IERR

      DOUBLE PRECISION SUMZ(NDOBSVYP),S2(NDOBSVYP),SUM,SUMY(NDOBSVZP),
     &  s(ndobsvzp)
      DOUBLE PRECISION OBSVYF(NDOBSVYP),SUMP,DIA

      DOUBLE PRECISION R(NDOBSVZP),PHI(NDOBSVYP)
     &  ,FPHI(NDOBSVYP)
     &  ,SZ(NDOBSVZP),SY(2*NDOBSVYP)
     &  ,FZ(NDOBSVZP),FY(NDOBSVYP),DPHI,DR,X,Y

      DOUBLE PRECISION OBSVZF(NDOBSVZP)

      DATA ICAL/0/

      save ical

      IF (IPINCIRC.EQ.0) THEN

        IF (IF1DIM.NE.2) THEN

C--- INTEGRATION ALONG HORIZONTAL AXIS Z

          IIY=0

          DO IY=1,nobsvy

            IIY=IIY+1
            OBSVYF(IIY)=OBSVY(IY)

            IF(MOBSVZ.GT.1) THEN

              iiz=0
              DO IZ=1,nobsvz
                iiz=iiz+1
                IOBSV=(IY-1)*NOBSVZ+IZ
                OBSVZF(IIZ)=OBSVZ(IZ)
                S(iiz)=spec(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(kfreq-1)))
              ENDDO      !IZ

              CALL FSPLINDX_omp(OBSVDZ,S,NOBSVZ,0.0D0,0.0D0,S2)

              SUMZ(IIY)=0.0d0

              DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ-1

                IOBSV=(IY-1)*NOBSVZ+IZ

                DSUM=
     &            OBSVDZ*0.5D0*(s(iz)+s(iz+1))-
     &            OBSVDZ**3/24.D0*(s2(iz)+s2(iz+1))

                IF(
     &              (IWSOUR.NE.ISOUR.OR.IWFREQ.NE.kfreq)
     &              .AND.
     &              DSUM.LT.0.0) THEN
                  WRITE(LUNGFO,*)
                  WRITE(LUNGFO,*)
                  WRITE(LUNGFO,*)'*** WARNING SR blende_omp ***'
                  WRITE(LUNGFO,*)
     &              'SPLINE INTEGRATION FAILED, RESULTS NOT RELIABLE'
                  WRITE(LUNGFO,*)'SOURCE POINT AND PHOTON ENERGY:'
     &              ,ISOUR,SNGL(FREQ(kfreq))
                  WRITE(LUNGFO,*)
                  WRITE(LUNGFO,*)
                  IWSOUR=ISOUR
                  IWFREQ=kfreq
                  IW_BLENF=1
                ENDIF !IWSOUR

                SUMZ(IIY)=SUMZ(IIY)+DSUM

              ENDDO   !IZ

            ELSE !MOBSVZ.GT.1

              IOBSV=(IY-1)*NOBSVZ+(NOBSVZ-MOBSVZ)/2+1
              SUMZ(IIY)=OBSVDZ*spec(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(kfreq-1)))

            ENDIF

          ENDDO !IY

        ELSE  !IF1DIM.EQ.2

C--- INTEGRATION ALONG HORIZONTAL AXIS Z

          IIY=0
          DO IY=1,nobsvy

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

            SUMZ(IIY)=DIA*spec(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(kfreq-1)))

          ENDDO !IY

        ENDIF !IF1DIM.EQ.2

C--- INTEGRATION ALONG VERTICAL AXIS Y

        CALL FSPLINDX_omp(OBSVDY,SUMZ,nobsvy,0.0D0,0.0D0,S2)

        IF(MOBSVY.GT.1) THEN

          SUM=0.0d0

          DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY-1

            DSUM=
     &        OBSVDY*0.5D0
     &        *(SUMZ(IY)+SUMZ(IY+1))
     &        -OBSVDY**3/24.D0
     &        *(S2(IY)+S2(IY+1))

            IF(
     &          (IWSOUR.NE.ISOUR.OR.IWFREQ.NE.kfreq)
     &          .AND.
     &          DSUM.LT.0.0) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** WARNING SR blende_omp ***'
              WRITE(LUNGFO,*)
     &          'SPLINE INTEGRATION FAILED, RESULTS NOT RELIABLE'
              WRITE(LUNGFO,*)'SOURCE POINT AND PHOTON ENERGY:'
     &          ,ISOUR,SNGL(FREQ(kfreq))
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              IWSOUR=ISOUR
              IWFREQ=kfreq
              IW_BLENF=1
            ENDIF !IWSOUR

            SUM=SUM+DSUM

          ENDDO

        ELSE IF (IF1DIM.EQ.2) THEN

          SUM=PI1*PINR*SUMZ(MMEDGEY+1)/2.0D0

        ELSE

          SUM=OBSVDY*SUMZ(MMEDGEY+1)

        ENDIF

      ELSE  !IPINCIRC

        IF (IRPHI.NE.0) THEN !INTEGRATION WITH RESPECT TO POLAR COORDINATES
          print*,"*** Obsolete option IRPHI not implemented for OMP version of WAVE ***"
          write(lungfo,*)"*** Obsolete option IRPHI not implemented for OMP version of WAVE ***"
          stop
        ELSE  !IRPHI

          DO IOBSV=1,NOBSV
            fphir_th(IOBSV)=spec(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(kfreq-1)))
          ENDDO !IOBSV

          CALL CIRCPIN_omp(NOBSVZ,NOBSVY,MOBSVZ,MOBSVY,SUM,SUMP,
     &      isour,kfreq,IERR)

          IF (IERR.NE.0) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'SR CIRCPIN_omp HAS BEEN CALLED BY SR blende_omp with error'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'SR CIRCPIN_omp HAS BEEN CALLED BY SR blende_omp with error'
            WRITE(6,*)
          ENDIF

        ENDIF !IRPHI

      ENDIF !PINCIRC

      wflux(ISOUR+NSOURCE*(kfreq-1))=SUM

      DO IY=1,NOBSVY
        DO IZ=1,NOBSVZ
          IOBSV=(IY-1)*NOBSVZ+IZ
          IOBFR=IOBSV+NOBSV*(kfreq-1)
          spectot(IOBFR)=spectot(IOBFR)+
     &      spec(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(kfreq-1)))
        ENDDO   !IZ
      ENDDO   !IY

      wfluxt(kfreq)=wfluxt(kfreq)+wflux(ISOUR+NSOURCE*(kfreq-1))

      RETURN
      END
