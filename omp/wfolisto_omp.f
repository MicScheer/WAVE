*CMZ :          03/10/2024  14.46.43  by  Michael Scheer
*CMZ :  3.07/00 13/03/2019  13.41.23  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.51/02 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.34/09 26/09/2001  15.58.40  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  16.27.20  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.15/00 08/05/2000  17.44.33  by  Michael Scheer
*CMZ :  2.13/10 14/04/2000  17.21.59  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.36  by  Michael Scheer
*CMZ :  2.13/03 12/01/2000  16.31.33  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.43.03  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.22.18  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  18.03.39  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.11  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE wfolisto_omp(ISTOK,kfreq)
*KEEP,gplhint.
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

      use wobsvmod

C--- FOLD SPECTRAL DENSITY IN PINHOLE

C    THE INTENSITY IS GIVEN ANALYTICALLY BY THE COEFFICINETS COFOLD
C    THE KERNEL IS GIVEN BY THE COEFFICIENTS GCOEFH, GCOEFV OF THE FOURIER
C    EXPANSION

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

      INTEGER IZA,IZE,IYA,IYE,kfreq,ISOUR,IY,IZ,IMASH,INZ
     &       ,INY,ISY,ISZ,I,J,M,KNZ,KNY,K,MARG
     &       ,ISMASH,IMASHZ,IMASHY,II,ISTOK,IGZY,NF,NFOLD,IFAIL


      DOUBLE PRECISION SUM,ARG,ZKZ0,ZKDZ,YKY0,YKDY,XKM1,DZ,DY

      DOUBLE PRECISION GSNZ(NGCOEFP,ngfourz),GCSZ(NGCOEFP,ngfourz)
      DOUBLE PRECISION GSNY(NGCOEFP,ngfoury),GCSY(NGCOEFP,ngfoury)

c      DOUBLE PRECISION ZINTK(NGCOEFP,4,NDMASHZP),ZINTKS(NGCOEFP,NDMASHZP)
c      DOUBLE PRECISION YINTK(NGCOEFP,4,NDMASHYP),YINTKS(NGCOEFP,NDMASHYP)

      DOUBLE PRECISION ZINTK(max(ngfourz,ngfoury),4,nobsvz),ZINTKS(max(ngfourz,ngfoury),nobsvz)
      DOUBLE PRECISION YINTK(max(ngfourz,ngfoury),4,NDMASHYP),YINTKS(max(ngfourz,ngfoury),NDMASHYP)

      DOUBLE PRECISION DCPZ,DCZ,DZXKM1,DZXKM12,GCPZ,GCZ,GSZ,GSPZ
      DOUBLE PRECISION DCPY,DCY,DYXKM1,DYXKM12,GCPY,GCY,GSY,GSPY

      DOUBLE PRECISION GZ(NGCOEFP*LIDIMP),GY(NGCOEFP*LIDIMP)
      EQUIVALENCE (GZ,GCOEFH),(GY,GCOEFV)

      ISOUR=ISIGSTO

      IF (IFOLD.EQ.1) THEN

        IYA=1+(NOBSVY-MOBSVY)/2
        IYE=NOBSVY-(NOBSVY-MOBSVY)/2

        IF (NOBSVZ.GT.1) THEN

          IZA=1+(NOBSVZ-MOBSVZ)/2
          IZE=NOBSVZ-(NOBSVZ-MOBSVZ)/2

          NF=NOBSV*(kfreq-1)

          DO IY=1,NOBSVY

            DO IZ=1,NOBSVZ
              IMASH=IZ+(IY-1)*NOBSVZ
c              WOBS1(IZ)=STOKES(ISTOK,IMASH+NF)
              x_th(iz)=obsvz(iz)
              WOBSv1_th(IZ)=STOKES(ISTOK,IMASH+NF)
c              if (iz.eq.nobsvz/2+1.and.iy.eq.nobsvy/2+1) then
c                print*,"iz,iy,imash+nf:",iz,iy,imash
c                print*,"z, stokes:",x_th(iz),wobsv1_th(iz)
c              endif
            ENDDO !IZ

            CALL UTIL_FOLD_FUNCTION_GAUSS_omp(
     &        NOBSVZ,WSIGZ(ISOUR),DGSIGZ(ISOUR))

c            CALL UTIL_FOLD_FUNCTION_GAUSS_omp(
c     &        NOBSVZ,OBSVZ,WOBS1,WSIGZ(ISOUR),DGSIGZ(ISOUR),WOBS2,
c     &        WOBS3,WOBS4,WOBS5,WOBS6,WOBS7)

            DO IZ=IZA,IZE
              IMASH=IZ+(IY-1)*NOBSVZ
c              STOKESF(ISTOK,IMASH+NF)=WOBS2(IZ)
              STOKESF(ISTOK,IMASH+NF)=WOBSv2_th(IZ)
c              if (iz.eq.nobsvz/2+1.and.iy.eq.nobsvy/2+1) then
c                print*,"iz,iy,imash+nf:",iz,iy,imash+nf
c                print*,"stokesf:",STOKESF(ISTOK,IMASH+NF)
c              endif
            ENDDO
          ENDDO   !IY

        ELSE !(NOBSVZ.GT.1)

          NF=NOBSV*(kfreq-1)

          DO IY=1,NOBSVY
            STOKESF(ISTOK,IY+NF)=STOKES(ISTOK,IY+NF)
          ENDDO

        ENDIF !(NOBSVZ.GT.1)

        DO IZ=1,NOBSVZ

          DO IY=1,NOBSVY
            IMASH=IZ+(IY-1)*NOBSVZ
c            WOBS1(IY)=STOKESF(ISTOK,IMASH+NF)
            x_th(iy)=obsvy(iy)
            WOBSv1_th(IY)=STOKESF(ISTOK,IMASH+NF)
c            if (iz.eq.nobsvz/2+1.and.iy.eq.nobsvy/2+1) then
c              print*,istok,iz,iy,imash+nf,x_th(iy),wobsv1_th(iy)
c            endif
          ENDDO   !IY

          CALL UTIL_FOLD_FUNCTION_GAUSS_omp(
     &      NOBSVY,WSIGy(ISOUR),DGSIGy(ISOUR))

c          CALL UTIL_FOLD_FUNCTION_GAUSS_omp(
c     &      NOBSVY,OBSVY,WOBS1,WSIGY(ISOUR),DGSIGY(ISOUR),WOBS2,
c     &      WOBS3,WOBS4,WOBS5,WOBS6,WOBS7)

          DO IY=IYA,IYE
            IMASH=IZ+(IY-1)*NOBSVZ
c            STOKESF(ISTOK,IMASH+NF)=WOBS2(IY)
            STOKESF(ISTOK,IMASH+NF)=WOBSv2_th(IY)
          ENDDO

        ENDDO  !IZ

      ELSE IF (IFOLD.EQ.-2) THEN

C--- ESTABLISH EDGE AROUND PINHOLE

        MOBSVZ=MOBSVZ+2*MMEDGEZ
        MOBSVY=MOBSVY+2*MMEDGEY
        MOBSV=MOBSVZ*MOBSVY

        IZA=1+(NOBSVZ-MOBSVZ)/2
        IZE=NOBSVZ-(NOBSVZ-MOBSVZ)/2

        IYA=1+(NOBSVY-MOBSVY)/2
        IYE=NOBSVY-(NOBSVY-MOBSVY)/2

        DZ=OBSVDZ
        DY=OBSVDY

C08052000      KNZ=DSIGZ(ISOUR)/DZ    !NUMBER OF ADJACENT MASHES
C08052000      KNY=DSIGY(ISOUR)/DY

        KNZ=IZA-1
        KNY=IYA-1

        KNZ=KNZ*2   !BOTH SIDES OF MASH POINT
        KNY=KNY*2

        IF (IF1DIM.EQ.1) KNZ=-99

        IF (KNZ.GT.NDMASHZ) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*) '*** ERROR IN wfolisto_omp ***'
          WRITE(LUNGFO,*) 'DIMENSION EXCEEDED. INCREASE PARAMETER'
          WRITE(LUNGFO,*) 'NDMASHZP IN CMPARA.CMN'
          WRITE(6,*) '*** ERROR IN wfolisto_omp ***'
          STOP
        ENDIF

        IF (KNY.GT.NDMASHY) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*) '*** ERROR IN wfolisto_omp ***'
          WRITE(LUNGFO,*) 'DIMENSION EXCEEDED. INCREASE PARAMETER'
          WRITE(LUNGFO,*) 'NDMASHYP IN CMPARA.CMN'
          WRITE(6,*) '*** ERROR IN wfolisto_omp ***'
          STOP
        ENDIF

C-- SET UP SIN AND COS ARRAYS

        ZKZ0=XKGAUSS(1,ISOUR)
        ZKDZ=XKGAUSS(1,ISOUR)*OBSVDZ

        DO INZ=1,KNZ+1
          DO M=1,NGFOURZ

C--- ATTENTION: IT'S NOT CLEAR WHAT'S RIGHT, BUT DOES IT MATTER
C             IF FOLDING FUNCTION IS SYMMETRIC ??

C        MARG=(M-1)*ABS(INZ-KNZ/2-1)
            MARG=(M-1)*(INZ-KNZ/2-1)

            ARG=DFLOAT(MARG)*ZKDZ
            GSNZ(M,INZ)=GCOEFH(M,ISOUR)*DSIN(ARG)
            GCSZ(M,INZ)=GCOEFH(M,ISOUR)*DCOS(ARG)

          ENDDO !M
        ENDDO !INZ

        YKY0=YKGAUSS(1,ISOUR)
        YKDY=YKGAUSS(1,ISOUR)*OBSVDY

        DO INY=1,KNY+1
          DO M=1,NGFOURY

C--- ATTENTION: IT'S NOT CLEAR WHAT'S RIGHT, BUT DOES IT MATTER
C             IF FOLDING FUNCTION IS SYMMETRIC ??

C        MARG=(M-1)*ABS(INY-KNY/2-1)
            MARG=(M-1)*(INY-KNY/2-1)

            ARG=DFLOAT(MARG)*YKDY
            GSNY(M,INY)=GCOEFV(M,ISOUR)*DSIN(ARG)
            GCSY(M,INY)=GCOEFV(M,ISOUR)*DCOS(ARG)

          ENDDO !M
        ENDDO !INY

C--- CALCULATE ELEMENTARY INTEGRALS

        DO INZ=1,KNZ  !LOOP OVER ADJACENT INTERVALLS
          DO M=1,NGFOURZ

            XKM1=DFLOAT(M-1)*ZKZ0

            GCPZ=GCSZ(M,INZ+1)
            GCZ=GCSZ(M,INZ)
            GSPZ=GSNZ(M,INZ+1)
            GSZ=GSNZ(M,INZ)

            DCPZ=GCPZ*DZ
            DCZ=GCZ*DZ
            DZXKM1=DZ*XKM1
            DZXKM12=DZXKM1*DZXKM1

            IF(M.EQ.1) THEN
              ZINTK(M,1,INZ)=DCPZ
              ZINTK(M,2,INZ)=DCPZ/2.
              ZINTK(M,3,INZ)=DCPZ/3.
              ZINTK(M,4,INZ)=DCPZ/4.

            ELSE

              ZINTK(M,1,INZ)=(GSPZ-GSZ)/XKM1
              ZINTK(M,2,INZ)=(GSPZ*DZXKM1+GCPZ-GCZ)/(DZXKM1*XKM1)
              ZINTK(M,3,INZ)=
     &          ((DZXKM12-2.)*GSPZ+2.*GSZ+2.*GCPZ*DZXKM1)
     &          /(DZXKM12*XKM1)
              ZINTK(M,4,INZ)=
     &          (3.*(DZXKM12-2.)*GCPZ+(DZXKM12-6.)*GSPZ*DZXKM1+6.
     &          *GCZ)/(DZXKM12*DZXKM1*XKM1)

            ENDIF
          ENDDO   !M

        ENDDO !INZ

        DO INY=1,KNY
          DO M=1,NGFOURY

            XKM1=DFLOAT(M-1)*YKY0

            GCPY=GCSY(M,INY+1)
            GCY=GCSY(M,INY)
            GSPY=GSNY(M,INY+1)
            GSY=GSNY(M,INY)

            DCPY=GCPY*DY
            DCY=GCY*DY
            DYXKM1=DY*XKM1
            DYXKM12=DYXKM1*DYXKM1

            IF(M.EQ.1) THEN
              YINTK(M,1,INY)=DCPY
              YINTK(M,2,INY)=DCPY/2.
              YINTK(M,3,INY)=DCPY/3.
              YINTK(M,4,INY)=DCPY/4.

            ELSE

              YINTK(M,1,INY)=(GSPY-GSY)/XKM1
              YINTK(M,2,INY)=(GSPY*DYXKM1+GCPY-GCY)/(DYXKM1*XKM1)
              YINTK(M,3,INY)=
     &          ((DYXKM12-2.)*GSPY+2.*GSY+2.*GCPY*DYXKM1)
     &          /(DYXKM12*XKM1)
              YINTK(M,4,INY)=
     &          (3.*(DYXKM12-2.)*GCPY+(DYXKM12-6.)*GSPY*DYXKM1+6.
     &          *GCY)/(DYXKM12*DYXKM1*XKM1)

            ENDIF
          ENDDO   !M

        ENDDO !INY


C--- SUM UP ELEMENTARY INTEGRALS (I.E. INNER PART OF INTEGRATION)

        DO K=1,KNZ
          DO I=1,4
            ZINTKS(I,K)=0.0
            DO M=1,NGFOURZ

              ZINTKS(I,K)=ZINTKS(I,K)+ZINTK(M,I,K)

            ENDDO !M
          ENDDO !I
        ENDDO !K

        DO K=1,KNY
          DO I=1,4
            YINTKS(I,K)=0.0
            DO M=1,NGFOURY

              YINTKS(I,K)=YINTKS(I,K)+YINTK(M,I,K)

            ENDDO !M
          ENDDO !I
        ENDDO !K

        IF (IF1DIM.EQ.0) THEN
          II=4
        ELSE
          KNZ=1
          II=1
          ZINTKS(1,1)=1.
        ENDIF

C--- LOOP OVER ALL MASHES INSIDE PINHOLE

        DO IY=IYA,IYE
          DO IZ=IZA,IZE

            IMASH=IZ+(IY-1)*NOBSVZ

C- LOOP OVER ALL ADJACENT MASHES AND SUM UP THE CONTRIBUTIONS

            SUM=0.0

            DO ISY=1,KNY

              IMASHY=IY-KNY/2+ISY-1

              DO ISZ=1,KNZ

                IMASHZ=IZ-KNZ/2+ISZ-1
                ISMASH=IMASHZ+(IMASHY-1)*NOBSVZ

C--- LOOP OVER ALL INDICIES OF BICUBIC SPLINE

                DO I=1,II
                  DO J=1,4

                    SUM=SUM+
     &                COFOLD(I,J,ISMASH)*ZINTKS(I,ISZ)*YINTKS(J,ISY)

                  ENDDO !J
                ENDDO !I

              ENDDO !ISZ
            ENDDO !ISY

            STOKESF(ISTOK,IMASH+NOBSV*(kfreq-1))=SUM

          ENDDO !IZ
        ENDDO !IY

      ELSE IF (IFOLD.EQ.-1) THEN

        IYA=1+(NOBSVY-MOBSVY)/2
        IYE=NOBSVY-(NOBSVY-MOBSVY)/2

        IF (NOBSVZ.GT.1) THEN

          IZA=1+(NOBSVZ-MOBSVZ)/2
          IZE=NOBSVZ-(NOBSVZ-MOBSVZ)/2

          ZKZ0=XKGAUSS(1,ISOUR)
          NFOLD=NINT(DSIGZ(ISOUR)/OBSVDZ)

          NF=NOBSV*(kfreq-1)
          IGZY=1+NGCOEFP*(ISOUR-1)

          DO IY=1,NOBSVY

            DO IZ=1,NOBSVZ
              IMASH=IZ+(IY-1)*NOBSVZ
              WOBS1(IZ)=STOKES(ISTOK,IMASH+NF)
            ENDDO !IZ

            CALL UTIL_FOLD_FOURIER(OBSVZ,WOBS1,NOBSVZ,NFOLD
     &        ,GZ(IGZY),NGFOURZ,WOBS2,WOBS3,WOBS4,WOBS5,WOBS6,WOBS7,IFAIL)

            IF (IFAIL.NE.0) WRITE(LUNGFO,*)
     &        '*** WARNING IN wfolisto_omp: FAILURE IN UTIL_FOLD_FOURIER ***'

            DO IZ=IZA,IZE
              IMASH=IZ+(IY-1)*NOBSVZ
              STOKESF(ISTOK,IMASH+NF)=WOBS2(IZ)
            ENDDO

          ENDDO   !IY

        ELSE !(NOBSVZ.GT.1)

          NF=NOBSV*(kfreq-1)
          IGZY=1+NGCOEFP*(ISOUR-1)
          DO IY=1,NOBSVY
            STOKESF(ISTOK,IY+NF)=STOKES(ISTOK,IY+NF)
          ENDDO

        ENDIF !(NOBSVZ.GT.1)

        YKY0=YKGAUSS(1,ISOUR)
        NFOLD=NINT(DSIGY(ISOUR)/OBSVDY)
        DO IZ=1,NOBSVZ

          DO IY=1,NOBSVY
            IMASH=IZ+(IY-1)*NOBSVZ
            WOBS1(IY)=STOKESF(ISTOK,IMASH+NF)
          ENDDO   !IY

          CALL UTIL_FOLD_FOURIER(OBSVY,WOBS1,NOBSVY,NFOLD
     &      ,GY(IGZY),NGFOURY,WOBS2,WOBS3,WOBS4,WOBS5,WOBS6,WOBS7,IFAIL)

          IF (IFAIL.NE.0) WRITE(LUNGFO,*)
     &      '*** WARNING IN wfolisto_omp: FAILURE IN UTIL_FOLD_FOURIER ***'

          DO IY=IYA,IYE
            IMASH=IZ+(IY-1)*NOBSVZ
            STOKESF(ISTOK,IMASH+NF)=WOBS2(IY)
          ENDDO

        ENDDO  !IZ

      ENDIF !IFOLD

c      print*,"int stokesf(1,18):",stokesf(1,18)

      RETURN
      END
