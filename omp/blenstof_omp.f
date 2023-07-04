*CMZ :  3.08/01 04/04/2019  08.55.47  by  Michael Scheer
*CMZ :  3.07/00 16/03/2019  15.27.40  by  Michael Scheer
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
      SUBROUTINE BLENSTOF_omp(ISTOK,kfreq)
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

      use circpinmod

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

      INTEGER kfreq,ISTOK,IY,IZ,IOBSV,IIY,IIZ
     &  ,ICAL,IR,IP,MR,MP,KDUM,IERR

      DOUBLE PRECISION SUMZ(NDOBSVYP),S2(NDOBSVYP),SUM,SUMY(NDOBSVZP)
      DOUBLE PRECISION OBSVYF(NDOBSVYP),DIA

      DOUBLE PRECISION R(NDOBSVZP),PHI(NDOBSVYP)
     &  ,FPHI(NDOBSVYP),s(ndobsvzp)
     &  ,SZ(NDOBSVZP),SY(NDOBSVYP)
     &  ,FZ(NDOBSVZP),FY(NDOBSVYP),DPHI,DR,X,Y

      DOUBLE PRECISION OBSVZF(NDOBSVZP)

      DATA ICAL/0/

      save ical

C--- TAKE INNER EDGE OF PINHOLE INTO ACCOUNT, I.E. SET MOBSVZ,MOBVY,MOBSV
C    TO ORIGINAL VALUES. THEY HAVE BEEN OVERWRITTEN IN SR WFOLINT

      MOBSVZ=MOBSVZ-2*MMEDGEZ
      MOBSVY=MOBSVY-2*MMEDGEY
      MOBSV=MOBSVZ*MOBSVY

c      print*,kfreq,ipincirc,irphi,if1dim,ical

      IF (IPINCIRC.EQ.0) THEN

        IF (IF1DIM.NE.2) THEN

C--- INTEGRATION ALONG HORIZONTAL AXIS Z

          IIY=0
          DO IY=(NOBSVY-MOBSVY-2*MMEDGEY)/2+1,
     &        (NOBSVY-MOBSVY-2*MMEDGEY)/2+MOBSVY+2*MMEDGEY

            IIY=IIY+1
            OBSVYF(IIY)=OBSVY(IY)

            iiz=0
            DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ
              iiz=iiz+1
              IOBSV=(IY-1)*NOBSVZ+IZ
              OBSVZF(IIZ)=OBSVZ(IZ)
              S(iiz)=STOKESF(ISTOK,IOBSV+NOBSV*(kfreq-1))
            ENDDO      !IZ

            CALL FSPLINDX_omp(OBSVDZ,S,MOBSVZ,0.D0,0.D0,S2)

            SUMZ(IIY)=0.0

            IF(MOBSVZ.GT.1) THEN

              iiz=0
              DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ-1

                iiz=iiz+1
                IOBSV=(IY-1)*NOBSVZ+IZ
                IOBFR=IOBSV+NOBSV*(kfreq-1)

                DSUM=
     &            OBSVDZ*0.5D0*(S(iiz)+s(iiz+1))
     &            -OBSVDZ**3/24.D0*(s2(iiz)+s2(iiz+1))

                SUMZ(IIY)=SUMZ(IIY)+DSUM

              ENDDO   !IZ

            ELSE

              IOBSV=(IY-1)*NOBSVZ+(NOBSVZ-MOBSVZ)/2+1
              SUMZ(IIY)=OBSVDZ*STOKESF(ISTOK,IOBSV+NOBSV*(kfreq-1))

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

            SUMZ(IIY)=DIA*STOKESF(ISTOK,IOBSV+NOBSV*(kfreq-1))

          ENDDO !IY

        ENDIF   !IF1DIM.EQ.2

C--- INTEGRATION ALONG VERTICAL AXIS Y

        CALL FSPLINDX_omp(OBSVDY,SUMZ,MOBSVY+2*MMEDGEY,0.0D0,0.0D0,S2)

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
          print*,"*** Obsolete option IRPHI not implemented for OMP version of WAVE ***"
          write(lungfo,*)"*** Obsolete option IRPHI not implemented for OMP version of WAVE ***"
          stop
        ELSE  !IRPHI

          DO IOBSV=1,NOBSV
            FPHIR_th(IOBSV)=STOKESF(ISTOK,IOBSV+NOBSV*(kfreq-1))
          ENDDO !IOBSV

          CALL CIRCPIN_omp(NOBSVZ,NOBSVY,MOBSVZ,MOBSVY,SUM,SUMP,
     &      -istok,kfreq,IERR)


        ENDIF !IRPHI

      ENDIF !PINCIRC

      WSTOKESF(ISTOK,kfreq)=SUM

      RETURN
      END
