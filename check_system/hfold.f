*CMZ :  4.00/14 30/12/2021  15.41.22  by  Michael Scheer
*CMZ :  4.00/13 07/12/2021  14.43.50  by  Michael Scheer
*CMZ :  3.03/01 07/10/2015  14.45.17  by  Michael Scheer
*CMZ :  3.02/06 24/04/2015  16.46.36  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  09.24.47  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  10.38.43  by  Michael Scheer
*CMZ :  2.57/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.51/02 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.51/00 26/05/2004  15.48.35  by  Michael Scheer
*CMZ :  2.47/09 27/05/2003  15.50.51  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  16.27.20  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  2.13/03 15/12/99  15.12.32  by  Michael Scheer
*CMZ :  2.13/02 14/12/99  16.33.52  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.28  by  Michael Scheer
*CMZ : 00.02/00 11/12/96  14.52.22  by  Michael Scheer
*CMZ : 00.01/12 16/10/96  15.22.18  by  Michael Scheer
*CMZ : 00.01/06 13/02/95  14.10.51  by  Michael Scheer
*CMZ : 00.01/04 30/01/95  14.13.04  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  16.53.40  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.52.27  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.49  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE HFOLD
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
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

C--- STORE RESULTS OF SPECTRUM CALCULATION ON HISTOGRAM FILE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,colli.
      include 'colli.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER IOB,ISOUR
      INTEGER IFO,IFOUR,ICYCLE

      DOUBLE PRECISION X

      DOUBLE PRECISION XFOLD,YFOLD,DXFOLD

      DOUBLE PRECISION FOUFUNX

      IF (NSOURCE.GE.10) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING SR HFOLD ***'
        WRITE(LUNGFO,*)'TOO MANY SOURCE POINTS'
        WRITE(LUNGFO,*)'HISTOGRAM IDENTIFIER WILL OVERLAP'
        WRITE(LUNGFO,*)'BE CAREFUL OR USE CORRESPONDING NTUPLE!!'
        WRITE(LUNGFO,*)
        WRITE(6,*)'*** WARNING SR HFOLD ***'
      ENDIF

C--- SCAN OBSERVATION POINTS

      X=OBSV(1,1)
      DO IOB=1,NOBSV
        IF(OBSV(1,IOB).NE.X) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN HFOLD ***'
          WRITE(LUNGFO,*)
     &      'X-POSITION OF OBSERVATION POINTS NOT IDENTICAL'
          WRITE(6,*)'*** ERROR IN HFOLD ***'
          WRITE(6,*)
     &      'X-POSITION OF OBSERVATION POINTS NOT IDENTICAL'
          STOP
        ENDIF
      ENDDO !IOB

C--- FOLDING FUNCTION

      DO ISOUR=1,NSOURCE

        call hbook1m(IDFOLFNY+ISOUR*100,'VERTICAL FOLDING FUNCTION'
     &    ,1000,-SNGL(DSIGY(ISOUR)),SNGL(DSIGY(ISOUR)),VMX)

        IF(IF1DIM.EQ.0) THEN

          call hbook1m(IDFOLFNZ+ISOUR*100,'HORIZONTAL FOLDING FUNCTION'
     &      ,1000,-SNGL(DSIGZ(ISOUR)),SNGL(DSIGZ(ISOUR)),VMX)

          DXFOLD=2.*DSIGZ(ISOUR)/1000.
          XFOLD=-DSIGZ(ISOUR)-DXFOLD/2.

          DO IFO=1,1000

            XFOLD=XFOLD+DXFOLD

            IF (IUSEM.EQ.0.AND.IPIN.NE.2) THEN

              IF (IFOLD.EQ.-1.OR.IFOLD.EQ.-2) THEN
                YFOLD=GCOEFH(1,ISOUR)
                DO IFOUR=2,NGFOURZ
                  YFOLD=YFOLD+GCOEFH(IFOUR,ISOUR)
     &              *DCOS(XKGAUSS(IFOUR-1,1)*XFOLD)
                ENDDO !IFOUR
              ELSE
                YFOLD=1.0D0/SQRT(2.0D0*PI1)/WSIGZ(ISOUR)*
     &            EXP(-(XFOLD/WSIGZ(ISOUR))**2/2.0D0)
              ENDIF !IFOLD

            ELSE  !IUSEM
              YFOLD=FOUFUNX(XFOLD,WSIGZ(ISOUR))
            ENDIF !IUSEM

            CALL hfillm(IDFOLFNZ+100*ISOUR,SNGL(XFOLD),0.,YFOLD)

          ENDDO   !IFO

        ENDIF !IF1DIM

        DXFOLD=2.*DSIGY(ISOUR)/1000.
        XFOLD=-DSIGY(ISOUR)-DXFOLD/2.

        DO IFO=1,1000
          XFOLD=XFOLD+DXFOLD

          IF (IUSEM.EQ.0.AND.IPIN.NE.2) THEN

            IF (IFOLD.EQ.-1.OR.IFOLD.EQ.-2) THEN
              YFOLD=GCOEFV(1,ISOUR)
              DO IFOUR=2,NGFOURY
                YFOLD=YFOLD+GCOEFV(IFOUR,ISOUR)
     &            *DCOS(YKGAUSS(IFOUR-1,1)*XFOLD)
              ENDDO !IFOUR
            ELSE
              YFOLD=1.0D0/SQRT(2.0D0*PI1)/WSIGY(ISOUR)*
     &          EXP(-(XFOLD/WSIGY(ISOUR))**2/2.0D0)
            ENDIF !IFOLD

          ELSE  !IUSEM
            YFOLD=FOUFUNX(XFOLD,WSIGY(ISOUR))
          ENDIF !IUSEM

          CALL hfillm(IDFOLFNY+100*ISOUR,SNGL(XFOLD),0.,YFOLD)

        ENDDO   !IFO

        if (if1dim.eq.0) then
          CALL MHROUT(IDFOLFNZ+100*ISOUR,ICYCLE,' ')
          CALL hdeletm(IDFOLFNZ+100*ISOUR)
        endif

        CALL MHROUT(IDFOLFNY+100*ISOUR,ICYCLE,' ')
        CALL hdeletm(IDFOLFNY+100*ISOUR)

      ENDDO   !ISOUR

      RETURN
      END
