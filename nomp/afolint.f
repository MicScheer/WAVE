*CMZ :          30/09/2024  14.57.39  by  Michael Scheer
*CMZ :  4.01/05 19/04/2024  12.22.35  by  Michael Scheer
*CMZ :  4.01/04 14/11/2023  13.46.13  by  Michael Scheer
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
*KEEP,GPLHINT.
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

C--- FOLD FIELD AMPLITUDE DENSITY IN PINHOLE

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

      DOUBLE PRECISION ZKZ0,YKY0,GZ(NGCOEFP*LIDIMP),GY(NGCOEFP*LIDIMP),rea(nobsv*nfreq)

c14.11.2023      if (icomp.eq.1) return

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
              rea(iobfr)=wobs2(iz)
c              if (icomp.eq.2.or.icomp.eq.3) then
c                reaima(icomp+2,ireim,iobfr)=WOBS2(IZ)
c              else if (icomp.eq.1) then
c                reaima(11,ireim,iobfr)=WOBS2(IZ)
c              else if (icomp.gt.5) then
c                reaima(icomp+6,ireim,iobfr)=WOBS2(IZ)
c              endif
            ENDDO

          ENDDO   !IY

        ELSE !(NOBSVZ.GT.1)

          NF=NOBSV*(IFREQ-1)

          DO IY=1,NOBSVY
            IMASH=1+(IY-1)*NOBSVZ
            IOBFR=imash+NOBSV*(IFREQ-1)
              rea(iobfr)=reaima(icomp,ireim,iobfr)
c            if (icomp.eq.2.or.icomp.eq.3) then
c              reaima(icomp+2,ireim,iobfr)=reaima(icomp,ireim,iobfr)
c            else if (icomp.eq.1) then
c              reaima(11,ireim,iobfr)=reaima(icomp,ireim,iobfr)
c            else if (icomp.gt.5) then
c              reaima(icomp+6,ireim,iobfr)=reaima(icomp,ireim,iobfr)
c            endif
          ENDDO

        ENDIF !(NOBSVZ.GT.1)

        DO IZ=1,NOBSVZ

          DO IY=1,NOBSVY
            IMASH=IZ+(IY-1)*NOBSVZ
            IOBFR=imash+NOBSV*(IFREQ-1)
c            wobs1(iy)=reaima(icomp,ireim,iobfr)
            wobs1(iy)=rea(iobfr)
          ENDDO   !IY

            CALL UTIL_FOLD_FUNCTION_GAUSS(
     &        NOBSVY,OBSVY,WOBS1,WSIGY(ISOUR),DGSIGY(ISOUR),WOBS2,
     &        WOBS3,WOBS4,WOBS5,WOBS6,WOBS7)

          DO IY=IYA,IYE
            IMASH=IZ+(IY-1)*NOBSVZ
            IOBFR=imash+NOBSV*(IFREQ-1)
            if (icomp.eq.2.or.icomp.eq.3) then
              reaima(icomp+2,ireim,iobfr)=WOBS2(IY)
            else if (icomp.eq.1) then
              reaima(11,ireim,iobfr)=WOBS2(IY)
            else if (icomp.gt.5) then
              reaima(icomp+6,ireim,iobfr)=WOBS2(IY)
            endif
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
              rea(iobfr)=WOBS2(IZ)
c              if (icomp.eq.2.or.icomp.eq.3) then
c                reaima(icomp+2,ireim,iobfr)=WOBS2(IZ)
c              else if (icomp.eq.1) then
c                reaima(11,ireim,iobfr)=WOBS2(IZ)
c              else if (icomp.gt.5) then
c                reaima(icomp+6,ireim,iobfr)=WOBS2(IZ)
c              endif
            ENDDO

          ENDDO   !IY

        ELSE !(NOBSVZ.GT.1)

          NF=NOBSV*(IFREQ-1)
          IGZY=1+NGCOEFP*(ISOUR-1)
          DO IY=1,NOBSVY
            IMASH=1+(IY-1)*NOBSVZ
            IOBFR=imash+NOBSV*(IFREQ-1)
            if (icomp.eq.2.or.icomp.eq.3) then
              reaima(icomp+2,ireim,iobfr)=reaima(icomp,ireim,iobfr)
            else if (icomp.eq.1) then
              reaima(11,ireim,iobfr)=reaima(icomp,ireim,iobfr)
            else if (icomp.gt.5) then
              reaima(icomp+6,ireim,iobfr)=reaima(icomp,ireim,iobfr)
            endif
          ENDDO

        ENDIF !(NOBSVZ.GT.1)

        YKY0=YKGAUSS(1,ISOUR)
        NFOLD=NINT(DSIGY(ISOUR)/OBSVDY)

        DO IZ=1,NOBSVZ

          DO IY=1,NOBSVY
            IMASH=IZ+(IY-1)*NOBSVZ
            IOBFR=imash+NOBSV*(IFREQ-1)
c            wobs1(iy)=reaima(icomp,ireim,iobfr)
            wobs1(iy)=rea(iobfr)
          ENDDO   !IY

          CALL UTIL_FOLD_FOURIER(OBSVY,WOBS1,NOBSVY,NFOLD
     &      ,GY(IGZY),NGFOURY,WOBS2,WOBS3,WOBS4,WOBS5,WOBS6,WOBS7,IFAIL)

          IF (IFAIL.NE.0) WRITE(LUNGFO,*)
     &      '*** WARNING IN AFOLINT: FAILURE IN UTIL_FOLD_FOURIER ***'

          DO IY=IYA,IYE
            IMASH=IZ+(IY-1)*NOBSVZ
            IOBFR=imash+NOBSV*(IFREQ-1)
            if (icomp.eq.2.or.icomp.eq.3) then
              reaima(icomp+2,ireim,iobfr)=WOBS2(IY)
            else if (icomp.eq.1) then
              reaima(11,ireim,iobfr)=WOBS2(IY)
            else if (icomp.gt.5) then
              reaima(icomp+6,ireim,iobfr)=WOBS2(IY)
            endif
          ENDDO

        ENDDO  !IZ

      ENDIF !IFOLD

      RETURN
      END
