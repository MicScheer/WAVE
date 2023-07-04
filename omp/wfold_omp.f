*CMZ :  3.07/00 15/03/2019  15.04.58  by  Michael Scheer
*CMZ :  3.03/02 21/03/2016  16.06.27  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/03 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.51/02 08/10/2009  09.58.11  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.46  by  Michael Scheer
*CMZ :  2.13/10 14/04/2000  14.26.49  by  Michael Scheer
*CMZ :  2.13/04 21/01/2000  14.54.46  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.22.28  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.38  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.16  by  Michael Scheer
*-- Author : Michael Scheer
*KEEP,wobsvmod.
      include 'wobsvmod.cmn'
*KEEP,circpinmod.
      include 'circpinmod.cmn'
*KEND.
      SUBROUTINE WFOLD_omp
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

      use omp_lib
      use ompmod
      use wobsvmod
      use circpinmod

C--- CALCULATES FOLDING OF PINHOLE INTENSITIY WITH ELECTRON PHASE SPACE
C    DISTRIBUTIONS (GAUSSIAN DISTRIBUTION)

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEND.

      INTEGER ISOUR,kfreq,IZ,IY,IOBSV
      !some variables do not appear correctly in the debugger, so replace them
      integer ith,kfold_th,kpincirc_th,krphi_th,kspecanaf_th

C--- CALCULATE FOURIER-COEFFICIENTS OF GAUSSIAN

      IF (ISPECANAF.NE.0) THEN
        print*,"*** WARNING: ISPECANAF not yet tested in WFOLD_OMP ***"
      ENDIF !ISPECANAF

      IF (IFOLD.NE.1) THEN
        print*,"*** WARNING: IFOLD.NE.1 not yet tested in WFOLD_OMP ***"
      ENDIF !ISPECANAF

      CALL WGFOUR

      SPECTOTF=0.0d0
      WFLUXTF=0.0d0

      kfold_th=ifold
      kpincirc_th=ipincirc
      krphi_th=irphi
      kspecanaf_th=ispecanaf

      DO iSOUR=1,NSOURCE

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMP& SHARED(mthreads,spec,specf,spectotf,obsv,obsvz,obsvy,wsigz,dgsigz,
!$OMP& iw_circ,if1dim,ipin,ispecanaf,ifold,ipincirc,irphi,wfluxf,wfluxtf,pinr,pincen)
!$OMP& FIRSTPRIVATE(nfreq,nobsv,isour,nsource,nobsvz,nobsvy,kfold_th,
!$OMP& kpincirc_th,krphi_th,kspecanaf_th,mobsv,mobsvy,mobsvz)

        ith=OMP_GET_THREAD_NUM()+1

        allocate(
     &    x_th(max(nobsvy,nobsvz)),
     &    wobsv1_th(max(nobsvy,nobsvz)),wobsv2_th(max(nobsvy,nobsvz)),
     &    wobsv3_th(max(nobsvy,nobsvz)),wobsv4_th(max(nobsvy,nobsvz)),
     &    wobsv5_th(max(nobsvy,nobsvz)),wobsv6_th(max(nobsvy,nobsvz)),
     &    wobsv7_th(max(nobsvy,nobsvz)))

        if (ipincirc.ne.0) then
          allocate(fphir_th(nobsv))
        endif

        ifold=kfold_th
        ipincirc=kpincirc_th
        irphi=krphi_th
        ispecanaf=kspecanaf_th

!$OMP DO

        DO kfreq=1,NFREQ

C--- CALCULATE 2D POLYNOMIALS OF INTENSITY DISTRIBUTION

          IF (kfold_th.EQ.-2) CALL WPOLY2(ISOUR,kfreq)

C--- PERFORM FOLDING

          CALL WFOLINT_omp(ith,ISOUR,kfreq)

C--- CALCULATE INTEGRATED INTENSITY

          IF (kspecanaf_th.NE.0) THEN
            CALL SPECANAF
          ENDIF !ISPECANAF

          CALL BLENDEF_omp(ISOUR,kfreq)

        ENDDO !kfreq

!$OMP END DO

        deallocate(x_th,
     &    wobsv1_th,wobsv2_th,wobsv3_th,wobsv4_th,wobsv5_th,
     &    wobsv6_th,wobsv7_th)

        if (ipincirc.ne.0) then
          deallocate(fphir_th)
        endif

!$OMP END PARALLEL

      ENDDO !ISOUR

      DO kfreq=1,NFREQ

        DO IY=1,NOBSVY
          DO IZ=1,NOBSVZ

            IF (IPINCIRC.EQ.0) THEN

              IF (
     &            IY.LT.(NOBSVY-MOBSVY)/2+1
     &            .OR.IY.GT.(NOBSVY-MOBSVY)/2+MOBSVY
     &            .OR.IZ.LT.(NOBSVZ-MOBSVZ)/2+1
     &            .OR.IZ.GT.(NOBSVZ-MOBSVZ)/2+MOBSVZ
     &            ) THEN

                SPECTOTF((IY-1)*NOBSVZ+IZ+NOBSV*(kfreq-1))=0.0

              ENDIF

            ELSE  !IPINCIRC

              IF (
     &            (OBSVZ(IZ)-PINCEN(3))**2
     &            +(OBSVY(IY)-PINCEN(2))**2
     &            -PINR**2
     &            .GT.1.D-10
     &            ) THEN

                SPECTOTF((IY-1)*NOBSVZ+IZ+NOBSV*(kfreq-1))=0.0

              ENDIF

            ENDIF !IPINCIRC

          ENDDO !IZ
        ENDDO !IY
      ENDDO !kfreq

      RETURN
      END
