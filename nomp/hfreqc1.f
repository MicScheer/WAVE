*CMZ :  4.00/15 07/04/2022  22.40.53  by  Michael Scheer
*CMZ :  4.00/14 30/12/2021  15.41.22  by  Michael Scheer
*CMZ :  4.00/13 07/12/2021  18.47.10  by  Michael Scheer
*CMZ :  4.00/04 23/08/2019  15.47.38  by  Michael Scheer
*CMZ :  3.04/01 03/04/2018  14.26.06  by  Michael Scheer
*CMZ :  3.03/02 05/01/2016  16.07.19  by  Michael Scheer
*CMZ :  3.03/00 24/09/2015  13.19.47  by  Michael Scheer
*CMZ :  3.02/06 08/06/2015  13.32.31  by  Michael Scheer
*CMZ :  3.02/04 06/03/2015  15.42.52  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/11 20/02/2013  16.35.31  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  11.38.54  by  Michael Scheer
*CMZ :  2.68/02 31/05/2012  14.03.31  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  10.38.44  by  Michael Scheer
*CMZ :  2.66/15 09/11/2010  16.06.02  by  Michael Scheer
*CMZ :  2.66/14 09/11/2010  14.11.04  by  Michael Scheer
*CMZ :  2.66/09 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.66/07 17/12/2009  16.10.22  by  Michael Scheer
*CMZ :  2.66/04 17/12/2009  16.08.16  by  Michael Scheer
*CMZ :  2.66/03 12/11/2009  16.27.11  by  Michael Scheer
*CMZ :  2.65/01 08/10/2009  09.58.11  by  Michael Scheer
*CMZ :  2.63/05 14/09/2009  15.19.42  by  Michael Scheer
*CMZ :  2.61/03 27/03/2007  13.23.57  by  Michael Scheer
*CMZ :  2.53/03 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.53/02 25/01/2005  18.15.44  by  Michael Scheer
*CMZ :  2.52/13 09/12/2004  13.09.09  by  Michael Scheer
*CMZ :  2.52/00 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.48/04 17/03/2004  13.47.13  by  Michael Scheer
*CMZ :  2.40/02 14/03/2002  15.46.20  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.45  by  Michael Scheer
*CMZ :  2.16/07 01/09/2000  14.36.07  by  Michael Scheer
*CMZ :  2.16/04 28/06/2000  17.43.04  by  Michael Scheer
*CMZ :  2.16/01 15/06/2000  15.47.11  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  16.17.53  by  Michael Scheer
*CMZ :  2.13/04 24/01/2000  16.35.38  by  Michael Scheer
*CMZ :  2.13/03 10/01/2000  17.32.03  by  Michael Scheer
*CMZ :  2.13/00 02/12/99  13.23.58  by  Michael Scheer
*CMZ :  1.03/00 16/01/98  11.16.49  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.27  by  Michael Scheer
*CMZ : 00.01/06 14/02/95  10.59.52  by  Michael Scheer
*CMZ : 00.01/05 31/01/95  17.02.40  by  Michael Scheer
*CMZ : 00.01/04 30/01/95  13.06.10  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  16.43.29  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.52.08  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.39  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE HFREQC1
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
*KEND.

C--- HISTOGRAMS FOR SPECTRA OF SINGLE OBSERVATION POINTS OR PINHOLE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEND.

      INTEGER ID,IFREQ,ICEN
      INTEGER ICYCLE,MFREQ

      REAL*4 FLOW,FHIG,DF
      DOUBLE PRECISION WEIGHT

      DF=FREQ(2)-FREQ(1)
      FLOW=FREQ(1)-DF/2.
      FHIG=FREQ(NFREQ)+DF/2.

      if (ifreq2p.eq.1.or.freqlow.eq.freqhig) then
        DF=freqhig-freqlow
        FLOW=freqlow-DF/2.
        FHIG=freqlow+DF/2.
      else if (ifreq2p.eq.-1) then
        DF=freqhig-freqlow
        FLOW=freqlow
        FHIG=freqhig
      endif

      IF (IPIN.ne.0) THEN

        ICEN=ICBRILL

        ID=ICFREQ
        MFREQ=max(1,NINT((FHIG-FLOW)/DF))
        call hbook1m(ID,'SELECTED FLUX DENSITY x 1.E-6 ',
     &    MFREQ,FLOW,FHIG,VMX)
        DO IFREQ=1,NFREQ,IHFREQ
          CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.
     &      ,SPECTOT(ICEN+NOBSV*(IFREQ-1))*1.0d-6)
        ENDDO   !NFREQ
        CALL MHROUT(ID,ICYCLE,' ')
        CALL hdeletm(ID)

        IF (IFOLD.ne.0) THEN
          ID=ICFREQF
          MFREQ=max(1,NINT((FHIG-FLOW)/DF))
          call hbook1m(ID,'SELECTED FLUX DENSITY (FOLDED) x 1.E-6 ',
     &      MFREQ,FLOW,FHIG,VMX)
          DO IFREQ=1,NFREQ,IHFREQ
            CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.
     &        ,SPECTOTF(ICEN+NOBSV*(IFREQ-1))*1.0d-6)
          ENDDO   !NFREQ
          CALL MHROUT(ID,ICYCLE,' ')
          CALL hdeletm(ID)
        ENDIF   !IFOLD

        IF (ISTOKES.NE.0) THEN

          ID=ICFRS0
          call hbook1m(ID,'SELECTED S0 x 1.E-6 '
     &      ,NFREQ,FLOW,FHIG,VMX)
          DO IFREQ=1,NFREQ,IHFREQ
            CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &        dble(STOKEC(1,IFREQ)*1.0d-6))
          ENDDO   !NFREQ
          CALL MHROUT(ID,ICYCLE,' ')
          CALL hdeletm(ID)

          ID=ICFRS1
          call hbook1m(ID,'SELECTED S1 x 1.E-6 '
     &      ,NFREQ,FLOW,FHIG,VMX)
          DO IFREQ=1,NFREQ,IHFREQ
            CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &        dble(STOKEC(2,IFREQ)*1.0d-6))
          ENDDO   !NFREQ
          CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

          ID=ICFRS2
          call hbook1m(ID,'SELECTED S2 x 1.E-6 '
     &      ,NFREQ,FLOW,FHIG,VMX)
          DO IFREQ=1,NFREQ,IHFREQ
            CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &        dble(STOKEC(3,IFREQ)*1.0d-6))
          ENDDO   !NFREQ
          CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

          ID=ICFRS3
          call hbook1m(ID,'SELECTED S3 x 1.E-6 '
     &      ,NFREQ,FLOW,FHIG,VMX)
          DO IFREQ=1,NFREQ,IHFREQ
            CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &        dble(stokEC(4,IFREQ)*1.E-6))
          ENDDO   !NFREQ
          CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

          ID=ICFRP
          call hbook1m(ID,'SELECTED P'
     &      ,NFREQ,FLOW,FHIG,VMX)
          DO IFREQ=1,NFREQ,IHFREQ
            WEIGHT=0.0d0
            IF (STOKEC(1,IFREQ).NE.0.0)
     &        WEIGHT=SQRT
     &        ((STOKEC(2,IFREQ)/STOKEC(1,IFREQ))**2
     &        +(STOKEC(3,IFREQ)/STOKEC(1,IFREQ))**2
     &        +(STOKEC(4,IFREQ)/STOKEC(1,IFREQ))**2)
            CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
          ENDDO   !NFREQ
          CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

          ID=ICFRP1
          call hbook1m(ID,'SELECTED P1'
     &      ,NFREQ,FLOW,FHIG,VMX)
          DO IFREQ=1,NFREQ,IHFREQ
            WEIGHT=0.0
            IF (STOKEC(1,IFREQ).NE.0.0)
     &        WEIGHT=STOKEC(2,IFREQ)/STOKEC(1,IFREQ)
            CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
          ENDDO   !NFREQ
          CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

          ID=ICFRP2
          call hbook1m(ID,'SELECTED P2'
     &      ,NFREQ,FLOW,FHIG,VMX)
          DO IFREQ=1,NFREQ,IHFREQ
            WEIGHT=0.0
            IF (STOKEC(1,IFREQ).NE.0.0)
     &        WEIGHT=STOKEC(3,IFREQ)/STOKEC(1,IFREQ)
            CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
          ENDDO   !NFREQ
          CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

          ID=ICFRP3
          call hbook1m(ID,'SELECTED P3'
     &      ,NFREQ,FLOW,FHIG,VMX)
          DO IFREQ=1,NFREQ,IHFREQ
            WEIGHT=0.0
            IF (STOKEC(1,IFREQ).NE.0.0)
     &        WEIGHT=STOKEC(4,IFREQ)/stokec(1,IFREQ)
            CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
          ENDDO   !NFREQ
          CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

          ID=ICFRP23
          call hbook1m(ID,'SELECTED P23'
     &      ,NFREQ,FLOW,FHIG,VMX)
          DO IFREQ=1,NFREQ,IHFREQ
            WEIGHT=0.0
            IF (stokec(1,IFREQ).NE.0.0)
     &        WEIGHT=SQRT((stokec(3,IFREQ)/stokec(1,IFREQ))**2
     &        +(stokec(4,IFREQ)/stokec(1,IFREQ))**2)
            CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
          ENDDO   !NFREQ
          CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRG3
              call hbook1m(ID,'SELECTED G3 x 1.e-6'
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                WEIGHT=0.0
                IF (stokec(1,IFREQ).NE.0.0)
     &            WEIGHT=stokec(4,IFREQ)/stokec(1,IFREQ)*stokec(4,IFREQ)*1.0d-6
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRG23
              call hbook1m(ID,'SELECTED G23 x 1.e-6'
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                WEIGHT=0.0
                IF (stokec(1,IFREQ).NE.0.0)
     &            WEIGHT=(stokec(3,IFREQ)/stokec(1,IFREQ)*stokec(3,IFREQ)
     &            +stokec(4,IFREQ)/stokec(1,IFREQ)*stokec(4,IFREQ))*1.0d-6
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              IF (IEFOLD.NE.0) THEN

                ID=ICFRS0E
                call hbook1m(ID,'SELECTED S0_E x 1.E-6 '
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &              dble(stokece(1,IFREQ)*1.e-6))
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRS1E
                call hbook1m(ID,'SELECTED S1_E x 1.E-6 '
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &              dble(stokece(2,IFREQ)*1.e-6))
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRS2E
                call hbook1m(ID,'SELECTED S2_E x 1.E-6 '
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &              dble(stokece(3,IFREQ)*1.e-6))
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRS3E
                call hbook1m(ID,'SELECTED S3_E x 1.E-6 '
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &              dble(stokece(4,IFREQ)*1.e-6))
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRPE
                call hbook1m(ID,'SELECTED P_E'
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (stokece(1,IFREQ).NE.0.0)
     &              WEIGHT=SQRT
     &              ((stokece(2,IFREQ)/stokece(1,IFREQ))**2
     &              +(stokece(3,IFREQ)/stokece(1,IFREQ))**2
     &              +(stokece(4,IFREQ)/stokece(1,IFREQ))**2)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRP1E
                call hbook1m(ID,'SELECTED P1_E'
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (stokece(1,IFREQ).NE.0.0)
     &              WEIGHT=stokece(2,IFREQ)/stokece(1,IFREQ)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRP2E
                call hbook1m(ID,'SELECTED P2_E'
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (stokece(1,IFREQ).NE.0.0)
     &              WEIGHT=stokece(3,IFREQ)/stokece(1,IFREQ)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRP3E
                call hbook1m(ID,'SELECTED P3_E'
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (stokece(1,IFREQ).NE.0.0)
     &              WEIGHT=stokece(4,IFREQ)/stokece(1,IFREQ)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRP23E
                call hbook1m(ID,'SELECTED P23_E'
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (stokece(1,IFREQ).NE.0.0)
     &              WEIGHT=SQRT((stokece(3,IFREQ)/stokece(1,IFREQ))**2
     &              +(stokece(4,IFREQ)/stokece(1,IFREQ))**2)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRG3E
                call hbook1m(ID,'SELECTED G3_E x 1.E-6 '
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (stokece(1,IFREQ).NE.0.0)
     &              WEIGHT=stokece(4,IFREQ)/stokece(1,IFREQ)*stokece(4,IFREQ)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight*1.0d-6)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRG23E
                call hbook1m(ID,'SELECTED G23_E x 1.E-6 '
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (stokece(1,IFREQ).NE.0.0)
     &              WEIGHT=stokece(3,IFREQ)/stokece(1,IFREQ)*stokece(3,IFREQ)
     &              +stokece(4,IFREQ)/stokece(1,IFREQ)*stokece(4,IFREQ)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight*1.0d-6)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

              ENDIF !IEFOLD

            IF (IFOLD.NE.0) THEN

              ID=ICFRS0F
              call hbook1m(ID,'SELECTED S0 (folded) x 1.E-6 '
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &            dble(stokECF(1,IFREQ)*1.e-6))
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRS1F
              call hbook1m(ID,'SELECTED S1 (folded) x 1.E-6 '
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &            dble(stokECF(2,IFREQ)*1.e-6))
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRS2F
              call hbook1m(ID,'SELECTED S2 (folded) x 1.E-6 '
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &            dble(stokECF(3,IFREQ)*1.e-6))
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRS3F
              call hbook1m(ID,'SELECTED S3 (folded) x 1.E-6 '
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &            dble(stokECF(4,IFREQ)*1.e-6))
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRPF
              call hbook1m(ID,'SELECTED P (folded)'
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                WEIGHT=0.0
                IF (STOKECF(1,IFREQ).NE.0.0)
     &            WEIGHT=SQRT
     &            ((STOKECF(2,IFREQ)/STOKECF(1,IFREQ))**2
     &            +(STOKECF(3,IFREQ)/STOKECF(1,IFREQ))**2
     &            +(STOKECF(4,IFREQ)/STOKECF(1,IFREQ))**2)
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRP1F
              call hbook1m(ID,'SELECTED P1 (folded)'
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                WEIGHT=0.0
                IF (STOKECF(1,IFREQ).NE.0.0)
     &            WEIGHT=STOKECF(2,IFREQ)/STOKECF(1,IFREQ)
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRP2F
              call hbook1m(ID,'SELECTED P2 (folded)'
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                WEIGHT=0.0
                IF (STOKECF(1,IFREQ).NE.0.0)
     &            WEIGHT=STOKECF(3,IFREQ)/STOKECF(1,IFREQ)
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRP3F
              call hbook1m(ID,'SELECTED P3 (folded)'
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                WEIGHT=0.0
                IF (STOKECF(1,IFREQ).NE.0.0)
     &            WEIGHT=STOKECF(4,IFREQ)/STOKECF(1,IFREQ)
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRP23F
              call hbook1m(ID,'SELECTED P23 (folded)'
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                WEIGHT=0.0
                IF (STOKECF(1,IFREQ).NE.0.0)
     &            WEIGHT=SQRT((STOKECF(3,IFREQ)/STOKECF(1,IFREQ))**2
     &            +(STOKECF(4,IFREQ)/STOKECF(1,IFREQ))**2)
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRG3F
              call hbook1m(ID,'SELECTED G3 (folded) x 1.e-6'
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                WEIGHT=0.0
                IF (STOKECF(1,IFREQ).NE.0.0)
     &            WEIGHT=STOKECF(4,IFREQ)/STOKECF(1,IFREQ)*STOKECF(4,IFREQ)*1.0d-6
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRG23F
              call hbook1m(ID,'SELECTED G23 (folded) x 1.e-6'
     &          ,NFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                WEIGHT=0.0
                IF (STOKECF(1,IFREQ).NE.0.0)
     &            WEIGHT=(STOKECF(3,IFREQ)/STOKECF(1,IFREQ)*STOKECF(3,IFREQ)
     &            +STOKECF(4,IFREQ)/STOKECF(1,IFREQ)*STOKECF(4,IFREQ))*1.0d-6
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              IF (IEFOLD.NE.0) THEN

                ID=ICFRS0EF
                call hbook1m(ID,'SELECTED S0_E (folded) x 1.E-6 '
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &              dble(stokECEF(1,IFREQ)*1.e-6))
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRS1EF
                call hbook1m(ID,'SELECTED S1_E (folded) x 1.E-6 '
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &              dble(stokECEF(2,IFREQ)*1.e-6))
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRS2EF
                call hbook1m(ID,'SELECTED S2_E (folded) x 1.E-6 '
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &              dble(stokECEF(3,IFREQ)*1.e-6))
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRS3EF
                call hbook1m(ID,'SELECTED S3_E (folded) x 1.E-6 '
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,
     &              dble(stokECEF(4,IFREQ)*1.e-6))
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRPEF
                call hbook1m(ID,'SELECTED P_E (folded)'
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (STOKECEF(1,IFREQ).NE.0.0)
     &              WEIGHT=SQRT
     &              ((STOKECEF(2,IFREQ)/STOKECEF(1,IFREQ))**2
     &              +(STOKECEF(3,IFREQ)/STOKECEF(1,IFREQ))**2
     &              +(STOKECEF(4,IFREQ)/STOKECEF(1,IFREQ))**2)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRP1EF
                call hbook1m(ID,'SELECTED P1_E (folded)'
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (STOKECEF(1,IFREQ).NE.0.0)
     &              WEIGHT=STOKECEF(2,IFREQ)/STOKECEF(1,IFREQ)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRP2EF
                call hbook1m(ID,'SELECTED P2_E (folded)'
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (STOKECEF(1,IFREQ).NE.0.0)
     &              WEIGHT=STOKECEF(3,IFREQ)/STOKECEF(1,IFREQ)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRP3EF
                call hbook1m(ID,'SELECTED P3_E (folded)'
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (STOKECEF(1,IFREQ).NE.0.0)
     &              WEIGHT=STOKECEF(4,IFREQ)/STOKECEF(1,IFREQ)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRP23EF
                call hbook1m(ID,'SELECTED P23_E (folded)'
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (STOKECEF(1,IFREQ).NE.0.0)
     &              WEIGHT=SQRT((STOKECEF(3,IFREQ)/STOKECEF(1,IFREQ))**2
     &              +(STOKECEF(4,IFREQ)/STOKECEF(1,IFREQ))**2)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRG3EF
                call hbook1m(ID,'SELECTED G3_E (folded) x 1.E-6 '
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (STOKECEF(1,IFREQ).NE.0.0)
     &              WEIGHT=STOKECEF(4,IFREQ)/STOKECEF(1,IFREQ)*STOKECEF(4,IFREQ)*1.0d-6
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
              CALL hdeletm(ID)

                ID=ICFRG23EF
                call hbook1m(ID,'SELECTED G23_E (folded) x 1.E-6 '
     &            ,NFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  WEIGHT=0.0
                  IF (STOKECEF(1,IFREQ).NE.0.0)
     &              WEIGHT=STOKECEF(3,IFREQ)/STOKECEF(1,IFREQ)*STOKECEF(3,IFREQ)
     &              +STOKECEF(4,IFREQ)/STOKECEF(1,IFREQ)*STOKECEF(4,IFREQ)
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,weight*1.0d-6)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
               call hdeletm(ID)

              ENDIF !IEFOLD

            ENDIF !IFOLD

          ENDIF !ISTOKES

          IF (IBRILL.NE.0) THEN

            ID=ICFRB0
            call hbook1m(ID,'SELECTED B0 x 1.E-12'
     &        ,MFREQ,FLOW,FHIG,VMX)
            DO IFREQ=1,NFREQ,IHFREQ
              CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillC(1,IFREQ))*1.E-12)
            ENDDO   !NFREQ
            CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

            ID=ICFRB1
            call hbook1m(ID,'SELECTED B1 x 1.E-12'
     &        ,MFREQ,FLOW,FHIG,VMX)
            DO IFREQ=1,NFREQ,IHFREQ
              CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillC(2,IFREQ))*1.E-12)
            ENDDO   !NFREQ
            CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

            ID=ICFRB2
            call hbook1m(ID,'SELECTED B2 x 1.E-12'
     &        ,MFREQ,FLOW,FHIG,VMX)
            DO IFREQ=1,NFREQ,IHFREQ
              CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillC(3,IFREQ))*1.E-12)
            ENDDO   !NFREQ
            CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

            ID=ICFRB3
            call hbook1m(ID,'SELECTED B3 x 1.E-12'
     &        ,MFREQ,FLOW,FHIG,VMX)
            DO IFREQ=1,NFREQ,IHFREQ
              CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillC(4,IFREQ))*1.E-12)
            ENDDO   !NFREQ
            CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

            IF (IEFOLD.NE.0) THEN

              ID=ICFRB0E
              call hbook1m(ID,'SELECTED B0_E x 1.E-12'
     &          ,MFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillCE(1,IFREQ))*1.E-12)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
               call hdeletm(ID)

              ID=ICFRB1E
              call hbook1m(ID,'SELECTED B1_E x 1.E-12'
     &          ,MFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillCE(2,IFREQ))*1.E-12)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
               call hdeletm(ID)

              ID=ICFRB2E
              call hbook1m(ID,'SELECTED B2_E x 1.E-12'
     &          ,MFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillCE(3,IFREQ))*1.E-12)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
               call hdeletm(ID)

              ID=ICFRB3E
              call hbook1m(ID,'SELECTED B3_E x 1.E-12'
     &          ,MFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillCE(4,IFREQ))*1.E-12)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
               call hdeletm(ID)

            ENDIF !IEFOLD

            IF (IFOLD.ne.2.and.ifold.ne.0) THEN

              ID=ICFRB0F
              call hbook1m(ID,'SELECTED B0 x 1.E-12 (FOLDED)'
     &          ,MFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillCF(1,IFREQ))*1.E-12)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRB1F
              call hbook1m(ID,'SELECTED B1 x 1.E-12 (FOLDED)'
     &          ,MFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillCF(2,IFREQ))*1.E-12)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRB2F
              call hbook1m(ID,'SELECTED B2 x 1.E-12 (FOLDED)'
     &          ,MFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillCF(3,IFREQ))*1.E-12)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              ID=ICFRB3F
              call hbook1m(ID,'SELECTED B3 x 1.E-12 (FOLDED)'
     &          ,MFREQ,FLOW,FHIG,VMX)
              DO IFREQ=1,NFREQ,IHFREQ
                CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillCF(4,IFREQ))*1.E-12)
              ENDDO   !NFREQ
              CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)

              IF (IEFOLD.NE.0) THEN

                ID=ICFRB0EF
                call hbook1m(ID,'SELECTED B0_E x 1.E-12 (FOLDED)'
     &            ,MFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillCEF(1,IFREQ))*1.E-12)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
               call hdeletm(ID)

                ID=ICFRB1EF
                call hbook1m(ID,'SELECTED B1_E x 1.E-12 (FOLDED)'
     &            ,MFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillCEF(2,IFREQ))*1.E-12)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
               call hdeletm(ID)

                ID=ICFRB2EF
                call hbook1m(ID,'SELECTED B2_E x 1.E-12 (FOLDED)'
     &            ,MFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillCEF(3,IFREQ))*1.E-12)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
               call hdeletm(ID)

                ID=ICFRB3EF
                call hbook1m(ID,'SELECTED B3_E x 1.E-12 (FOLDED)'
     &            ,MFREQ,FLOW,FHIG,VMX)
                DO IFREQ=1,NFREQ,IHFREQ
                  CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.,dble(brillCEF(4,IFREQ))*1.E-12)
                ENDDO   !NFREQ
                CALL MHROUT(ID,ICYCLE,' ')
               call hdeletm(ID)

              ENDIF !IEFOLD
            ENDIF !IFOLD.NE.2

          ENDIF !IBRILL

        ELSE    !IPIN

          ICEN=ICBRILL

          ID=ICFREQ
          MFREQ=max(1,NINT((FHIG-FLOW)/DF))
          call hbook1m(ID,'SELECTED FLUX DENSITY x 1.E-6 ',
     &      MFREQ,FLOW,FHIG,VMX)
          DO IFREQ=1,NFREQ,IHFREQ
            CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.
     &        ,specTOT(ICEN+NOBSV*(IFREQ-1))*1.0d-6)
          ENDDO   !NFREQ
          CALL MHROUT(ID,ICYCLE,' ')
          CALL hdeletm(ID)

          IF (IFOLD.NE.0) THEN
            ID=ICFREQF
            MFREQ=max(1,NINT((FHIG-FLOW)/DF))
            call hbook1m(ID,'SELECTED FLUX DENSITY  x 1.E-6 (FOLDED)',
     &        MFREQ,FLOW,FHIG,VMX)
            DO IFREQ=1,NFREQ,IHFREQ
              CALL hfillm(ID,SNGL(FREQ(IFREQ)),0.
     &          ,specTOTF(ICEN+NOBSV*(IFREQ-1))*1.0d-6)
            ENDDO   !NFREQ
            CALL MHROUT(ID,ICYCLE,' ')
            call hdeletm(ID)
          ENDIF   !IFOLD

        ENDIF    !IPIN

      RETURN
      END
