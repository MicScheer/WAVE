*CMZ :  3.07/00 08/03/2019  18.44.03  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.51/02 29/04/2010  11.46.31  by  Michael Scheer
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
      SUBROUTINE WFOLINT_omp(ith,ISOUR,kfreq)
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

      use wobsvmod

C--- FOLD SPECTRAL DENSITY IN PINHOLE

C    IFOLD.EQ.-2:
C    THE INTENSITY IS GIVEN ANALYTICALLY BY THE COEFFICINETS COFOLD
C    THE KERNEL IS GIVEN BY THE COEFFICIENTS GCOEFH, GCOEFV OF THE FOURIER
C    EXPANSION
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

      INTEGER IZA,IZE,IYA,IYE,kfreq,ISOUR,IY,IZ,IMASH,INZ
     &  ,INY,ISY,ISZ,I,J,M,KNZ,KNY,K,MARG
     &  ,ISMASH,IMASHZ,IMASHY,II,NF,NFOLD,IFAIL,IGZY

      integer ith,ifold_th

      DOUBLE PRECISION SUM,ARG,ZKZ0,ZKDZ,YKY0,YKDY,XKM1,DZ,DY,PI,YLY0

      DOUBLE PRECISION GSNZ(NGCOEFP,NDMASHZP),GCSZ(NGCOEFP,NDMASHZP)
      DOUBLE PRECISION GZ(NGCOEFP*LIDIMP),GY(NGCOEFP*LIDIMP)
      DOUBLE PRECISION GSNY(NGCOEFP,NDMASHYP),GCSY(NGCOEFP,NDMASHYP)
      DOUBLE PRECISION ZINTK(NGCOEFP,4,NDMASHZP),ZINTKS(NGCOEFP,NDMASHZP)
      DOUBLE PRECISION YINTK(NGCOEFP,4,NDMASHYP),YINTKS(NGCOEFP,NDMASHYP)

      DOUBLE PRECISION DCPZ,DCZ,DZXKM1,DZXKM12,GCPZ,GCZ,GSZ,GSPZ
      DOUBLE PRECISION DCPY,DCY,DYXKM1,DYXKM12,GCPY,GCY,GSY,GSPY

      EQUIVALENCE (GZ,GCOEFH),(GY,GCOEFV)

      DATA PI/3.141592653589793D0/

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
              x_th(iz)=obsvz(iz)
              wobsv1_th(IZ)=spec(ISOUR+NSOURCE*(IMASH-1+NF))
            ENDDO !IZ

            CALL UTIL_FOLD_FUNCTION_GAUSS_omp(
     &        NOBSVZ,WSIGZ(ISOUR),DGSIGZ(ISOUR))

            DO IZ=IZA,IZE
              IMASH=IZ+(IY-1)*NOBSVZ
              SPECF(ISOUR+NSOURCE*(IMASH-1+NF))=wobsv2_th(IZ)
            ENDDO

          ENDDO   !IY

        ELSE !(NOBSVZ.GT.1)

          NF=NOBSV*(kfreq-1)

          DO IY=1,NOBSVY
            SPECF(ISOUR+NSOURCE*((IY-1)*NOBSVZ+NF))=
     &        SPEC(ISOUR+NSOURCE*((IY-1)*NOBSVZ+NF))
          ENDDO

        ENDIF !(NOBSVZ.GT.1)

        DO IZ=1,NOBSVZ

          DO IY=1,NOBSVY
            IMASH=IZ+(IY-1)*NOBSVZ
            x_th(iy)=obsvy(iy)
            wobsv1_th(IY)=SPECF(ISOUR+NSOURCE*(IMASH-1+NF))
          ENDDO   !IY

          CALL UTIL_FOLD_FUNCTION_GAUSS_omp(
     &      NOBSVY,WSIGY(ISOUR),DGSIGY(ISOUR))

          DO IY=IYA,IYE
            IMASH=IZ+(IY-1)*NOBSVZ
            SPECF(ISOUR+NSOURCE*(IMASH-1+NF))=wobsv2_th(IY)
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

C      KNZ=DSIGZ(ISOUR)/DZ    !NUMBER OF ADJACENT MASHES
C      KNY=DSIGY(ISOUR)/DY

        KNZ=IZA-1
        KNY=IYA-1

        KNZ=KNZ*2   !BOTH SIDES OF MASH POINT
        KNY=KNY*2

        IF (IF1DIM.NE.0) KNZ=-99

        IF (KNZ.GT.NDMASHZ) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*) '*** ERROR IN WFOLINT_OMP ***'
          WRITE(LUNGFO,*) 'DIMENSION EXCEEDED. INCREASE PARAMETER'
          WRITE(LUNGFO,*) 'NDMASHZP IN CMPARA.CMN'
          WRITE(6,*) '*** ERROR IN WFOLINT_OMP ***'
          STOP
        ENDIF

        IF (KNY.GT.NDMASHY) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*) '*** ERROR IN WFOLINT_OMP ***'
          WRITE(LUNGFO,*) 'DIMENSION EXCEEDED. INCREASE PARAMETER'
          WRITE(LUNGFO,*) 'NDMASHYP IN CMPARA.CMN'
          WRITE(6,*) '*** ERROR IN WFOLINT_OMP ***'
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

            SPECF(ISOUR+NSOURCE*(IMASH-1+NOBSV*(kfreq-1)))=SUM

          ENDDO !IZ
        ENDDO !IY

      ELSE IF(IFOLD.EQ.-1) THEN

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
              wobsv1_th(IZ)=SPEC(ISOUR+NSOURCE*(IMASH-1+NF))
            ENDDO !IZ

            CALL UTIL_FOLD_FOURIER(OBSVZ,wobsv1_th,NOBSVZ,NFOLD
     &        ,GZ(IGZY),NGFOURZ,wobsv2_th,wobsv3_th,wobsv4_th,wobsv5_th,wobsv6_th,wobsv7_th,IFAIL)

            IF (IFAIL.NE.0) WRITE(LUNGFO,*)
     &        '*** WARNING IN WFOLINT_OMP: FAILURE IN UTIL_FOLD_FOURIER ***'

            DO IZ=IZA,IZE
              IMASH=IZ+(IY-1)*NOBSVZ
              SPECF(ISOUR+NSOURCE*(IMASH-1+NF))=wobsv2_th(IZ)
            ENDDO

          ENDDO   !IY

        ELSE !(NOBSVZ.GT.1)

          NF=NOBSV*(kfreq-1)
          IGZY=1+NGCOEFP*(ISOUR-1)
          DO IY=1,NOBSVY
            SPECF(ISOUR+NSOURCE*((IY-1)*NOBSVZ+NF))=
     &        SPEC(ISOUR+NSOURCE*((IY-1)*NOBSVZ+NF))
          ENDDO

        ENDIF !(NOBSVZ.GT.1)

        YKY0=YKGAUSS(1,ISOUR)
        YLY0=2.D0*PI/YKY0
        NFOLD=NINT(DSIGY(ISOUR)/OBSVDY)

        DO IZ=1,NOBSVZ

          DO IY=1,NOBSVY
            IMASH=IZ+(IY-1)*NOBSVZ
            wobsv1_th(IY)=SPECF(ISOUR+NSOURCE*(IMASH-1+NF))
          ENDDO   !IY

          CALL UTIL_FOLD_FOURIER(OBSVY,wobsv1_th,NOBSVY,NFOLD
     &      ,GY(IGZY),NGFOURY,wobsv2_th,wobsv3_th,wobsv4_th,wobsv5_th,wobsv6_th,wobsv7_th,IFAIL)

          IF (IFAIL.NE.0) WRITE(LUNGFO,*)
     &      '*** WARNING IN WFOLINT_OMP: FAILURE IN UTIL_FOLD_FOURIER ***'

          DO IY=IYA,IYE
            IMASH=IZ+(IY-1)*NOBSVZ
            SPECF(ISOUR+NSOURCE*(IMASH-1+NF))=wobsv2_th(IY)
          ENDDO

        ENDDO  !IZ

      ENDIF !IFOLD

c      print*,"Ende von wfolint_omp",ith

      RETURN
      END
