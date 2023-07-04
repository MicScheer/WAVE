*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/03 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.51/02 08/10/2009  09.58.11  by  Michael Scheer
*CMZ :  2.34/09 24/09/2001  16.48.01  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  16.27.20  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.15/00 08/05/2000  13.33.58  by  Michael Scheer
*CMZ :  2.14/02 27/04/2000  17.54.22  by  Michael Scheer
*CMZ :  2.13/10 14/04/2000  17.10.49  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.36  by  Michael Scheer
*CMZ :  2.13/03 20/12/99  17.39.24  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.43.03  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.21.35  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  18.07.26  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.22  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE AFOLINT(icomp,ireim,IFREQ)
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
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

C--- FOLD FIELDAMPLITUDE DENSITY IN PINHOLE

C    IFOLD.EQ.-2: not available
C    IFOLD.EQ.-1: SR UTIL_FOLD_FOURIER IS USED
C    IFOLD.EQ.1: SR UTIL_FOLD_FUNCTION_GAUSS IS USED

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEND.

      INTEGER IZA,IZE,IYA,IYE,IFREQ,IY,IZ,IMASH
     &       ,NF,NFOLD,IFAIL,IGZY,icomp,ireim,isour

      DOUBLE PRECISION ZKZ0,YKY0,GZ(NGCOEFP*LIDIMP),GY(NGCOEFP*LIDIMP)

      if (icomp.eq.1) return

      isour=nsource/2+1

      IF (IFOLD.EQ.1) THEN

        IYA=1+(NOBSVY-MOBSVY)/2
        IYE=NOBSVY-(NOBSVY-MOBSVY)/2

        IF (NOBSVZ.GT.1) THEN

          IZA=1+(NOBSVZ-MOBSVZ)/2
          IZE=NOBSVZ-(NOBSVZ-MOBSVZ)/2

          NF=NOBSV*(IFREQ-1)

          DO IY=1,NOBSVY

            DO IZ=1,NOBSVZ
              IMASH=IZ+(IY-1)*NOBSVZ
              IOBFR=imash+NOBSV*(IFREQ-1)
              wobs1(iz)=reaima(icomp,ireim,iobfr)
            ENDDO !IZ

            CALL UTIL_FOLD_FUNCTION_GAUSS(
     &        NOBSVZ,OBSVZ,WOBS1,WSIGZ(ISOUR),DGSIGZ(ISOUR),WOBS2,
     &        WOBS3,WOBS4,WOBS5,WOBS6,WOBS7)

            DO IZ=IZA,IZE
              IMASH=IZ+(IY-1)*NOBSVZ
              IOBFR=imash+NOBSV*(IFREQ-1)
              reaima(icomp+2,ireim,iobfr)=WOBS2(IZ)
            ENDDO

          ENDDO   !IY

        ELSE !(NOBSVZ.GT.1)

          NF=NOBSV*(IFREQ-1)

          DO IY=1,NOBSVY
            IMASH=1+(IY-1)*NOBSVZ
            IOBFR=imash+NOBSV*(IFREQ-1)
            reaima(icomp+2,ireim,iobfr)=reaima(icomp,ireim,iobfr)
          ENDDO

        ENDIF !(NOBSVZ.GT.1)

        DO IZ=1,NOBSVZ

          DO IY=1,NOBSVY
            IMASH=IZ+(IY-1)*NOBSVZ
            IOBFR=imash+NOBSV*(IFREQ-1)
            wobs1(iy)=reaima(icomp,ireim,iobfr)
          ENDDO   !IY

            CALL UTIL_FOLD_FUNCTION_GAUSS(
     &        NOBSVY,OBSVY,WOBS1,WSIGY(ISOUR),DGSIGY(ISOUR),WOBS2,
     &        WOBS3,WOBS4,WOBS5,WOBS6,WOBS7)

          DO IY=IYA,IYE
            IMASH=IZ+(IY-1)*NOBSVZ
            IOBFR=imash+NOBSV*(IFREQ-1)
            reaima(icomp+2,ireim,iobfr)=WOBS2(IY)
          ENDDO

        ENDDO  !IZ

      ELSE IF(IFOLD.EQ.-1) THEN

        IYA=1+(NOBSVY-MOBSVY)/2
        IYE=NOBSVY-(NOBSVY-MOBSVY)/2

        IF (NOBSVZ.GT.1) THEN

          IZA=1+(NOBSVZ-MOBSVZ)/2
          IZE=NOBSVZ-(NOBSVZ-MOBSVZ)/2

          ZKZ0=XKGAUSS(1,ISOUR)
          NFOLD=NINT(DSIGZ(ISOUR)/OBSVDZ)

          NF=NOBSV*(IFREQ-1)
          IGZY=1+NGCOEFP*(ISOUR-1)

          DO IY=1,NOBSVY

            DO IZ=1,NOBSVZ
              IMASH=IZ+(IY-1)*NOBSVZ
              IOBFR=imash+NOBSV*(IFREQ-1)
              wobs1(iz)=reaima(icomp,ireim,iobfr)
            ENDDO !IZ

            CALL UTIL_FOLD_FOURIER(OBSVZ,WOBS1,NOBSVZ,NFOLD
     &        ,GZ(IGZY),NGFOURZ,WOBS2,WOBS3,WOBS4,WOBS5,WOBS6,WOBS7,IFAIL)

            IF (IFAIL.NE.0) WRITE(LUNGFO,*)
     &        '*** WARNING IN AFOLINT: FAILURE IN UTIL_FOLD_FOURIER ***'

            DO IZ=IZA,IZE
              IMASH=IZ+(IY-1)*NOBSVZ
              IOBFR=imash+NOBSV*(IFREQ-1)
              reaima(icomp+2,ireim,iobfr)=WOBS2(IZ)
            ENDDO

          ENDDO   !IY

        ELSE !(NOBSVZ.GT.1)

          NF=NOBSV*(IFREQ-1)
          IGZY=1+NGCOEFP*(ISOUR-1)
          DO IY=1,NOBSVY
            IMASH=1+(IY-1)*NOBSVZ
            IOBFR=imash+NOBSV*(IFREQ-1)
            reaima(icomp+2,ireim,iobfr)=reaima(icomp,ireim,iobfr)
          ENDDO

        ENDIF !(NOBSVZ.GT.1)

        YKY0=YKGAUSS(1,ISOUR)
        NFOLD=NINT(DSIGY(ISOUR)/OBSVDY)

        DO IZ=1,NOBSVZ

          DO IY=1,NOBSVY
            IMASH=IZ+(IY-1)*NOBSVZ
            IOBFR=imash+NOBSV*(IFREQ-1)
            wobs1(iy)=reaima(icomp,ireim,iobfr)
          ENDDO   !IY

          CALL UTIL_FOLD_FOURIER(OBSVY,WOBS1,NOBSVY,NFOLD
     &      ,GY(IGZY),NGFOURY,WOBS2,WOBS3,WOBS4,WOBS5,WOBS6,WOBS7,IFAIL)

          IF (IFAIL.NE.0) WRITE(LUNGFO,*)
     &      '*** WARNING IN AFOLINT: FAILURE IN UTIL_FOLD_FOURIER ***'

          DO IY=IYA,IYE
            IMASH=IZ+(IY-1)*NOBSVZ
            IOBFR=imash+NOBSV*(IFREQ-1)
            reaima(icomp+2,ireim,iobfr)=WOBS2(IY)
          ENDDO

        ENDDO  !IZ

      ENDIF !IFOLD

      RETURN
      END
