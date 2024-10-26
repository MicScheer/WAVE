*CMZ :  3.07/00 08/03/2019  20.03.45  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.63/03 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.57/03 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.57/00 22/11/2005  11.19.13  by  Michael Scheer
*CMZ :  2.34/00 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.33/09 10/05/2001  18.02.55  by  Michael Scheer
*CMZ :  2.16/08 24/10/2000  11.19.57  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.43.04  by  Michael Scheer
*CMZ : 00.02/04 25/02/97  17.37.15  by  Michael Scheer
*CMZ : 00.01/09 01/09/95  13.03.01  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.45.26  by  Michael Scheer
*CMZ : 00.00/06 29/04/94  21.43.32  by  Michael Scheer
*CMZ : 00.00/05 29/04/94  20.18.16  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.15  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE circpin1d_omp(NY,MY,SPEC,RESULT,SUMP)
*KEEP,gplhint.
*KEND.

*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

      use wobsvmod

C---  SUBROUTINE TO INTEGRATE SPEC OVER CIRCULAR PINHOLE FOR IF1DIM OR VERT.
C     DIRECTION FOR IF1DIM=0
C     NZ,NY ARE ARRAY DIMENSIONS
C     INTEGRATION IS DONE FOR MY POINTS

      IMPLICIT NONE
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER NY,IY,MY,IIY,IIY1,NB

      DOUBLE PRECISION RESULT,SUMP
     &  ,SPEC(NDOBSVYP)
c     &  ,UF(NDOBSVYP)
c     &  ,UY(NDOBSVYP)
c     &  ,W1(NDOBSVYP),W2(NDOBSVYP),W3(NDOBSVYP),W4(NDOBSVYP),COEFF(NDOBSVYP)
     &  ,ANS1,ANS2,ANS3,ANS,YLO,YHI,FLO,FHI,YLO2,YLO3,YHI2,YHI3,R4,PINR2,
     &  F2LO,F2HI,SQRRYH2,SQRRYL2,SQLO,SQHI,ASINLO,ASINHI

      IF (MY.LT.2) THEN
        WRITE(6,*) '*** ERROR IN circpin1d_omp: MY TOO LOW ***'
        STOP
      ENDIF

C--- LOOP OF VERTICAL DIRECTION

      RESULT=0.0D0
      PINR2=PINR**2
      R4=PINR2**2

      IF (IF1DIM.NE.0) THEN

        DO IY=1,NY
          x_th(IY)=OBSVY(IY)-PINCEN(2)
          wobsv1_th(IY)=SPEC(IY)
        ENDDO !IY

        NB=(NY-MY)/2
        IIY=NY

      ELSE

        IIY=0
        DO IY=(NY-MY)/2+1,(NY-MY)/2+MY
          IIY=IIY+1
          x_th(IIY)=OBSVY(IY)-PINCEN(2)
          wobsv1_th(IIY)=SPEC(IIY)
        ENDDO !IY

        NB=0

      ENDIF !IF1DIM

      CALL UTIL_SPLINE_COEF_omp(iiy,-9999.0d0,-9999.0d0)
c      CALL UTIL_SPLINE_COEF(UY,UF,IIY,-9999.0d0,-9999.0d0,COEFF,W1,W2,W3,W4)

      DO IIY=NB+2,NB+MY

        IIY1=IIY-1

        YLO=x_th(IIY1)
        YHI=x_th(IIY)
        FLO=wobsv1_th(IIY1)
        FHI=wobsv1_th(IIY)
        F2LO=wobsv3_th(IIY1)
        F2HI=wobsv3_th(IIY)

        YLO2=YLO**2
        YHI2=YHI**2
        YLO3=YLO2*YLO
        YHI3=YHI2*YHI

        SQLO=PINR2-YLO2
        SQHI=PINR2-YHI2

        IF (SQLO.GT.0.0D0) THEN
          SQRRYL2=SQRT(SQLO)
        ELSE
          SQRRYL2=0.0D0
        ENDIF

        IF (SQHI.GT.0.0D0) THEN
          SQRRYH2=SQRT(SQHI)
        ELSE
          SQRRYH2=0.0D0
        ENDIF

        SQLO=YLO/PINR
        SQHI=YHI/PINR

        IF (SQLO.LE.-1.0D0) THEN
          ASINLO=-PI1/2.D0
        ELSE
          ASINLO=ASIN(SQLO)
        ENDIF

        IF (SQHI.GE.1.0D0) THEN
          ASINHI=PI1/2.D0
        ELSE
          ASINHI=ASIN(SQHI)
        ENDIF

        ans3=-15.0d0*((4.0d0*(2.0d0*yhi-ylo)*ylo+3.0d0*PINR2)*f2lo*yhi-24.0d0*(fhi
     &    *ylo-flo*yhi)+(4.0d0*(yhi-2.0d0*ylo)*yhi-3.0d0*PINR2)*f2hi*ylo)*asinhi
     &    *PINR2

        ans2=
     &    (
     &    2.0d0*((20.0d0*yhi2-25.0d0*yhi*ylo+8.0d0*ylo2)*sqrryl2*ylo2
     &    +(7.0d0*yhi2-20.0d0*yhi*ylo+10.0d0*ylo2)*sqrryh2*yhi2-8.0d0*(
     &    sqrryh2-sqrryl2)*r4)+((80.0d0*yhi2+35.0d0*yhi*ylo-32.0d0*ylo2)
     &    *sqrryl2-(43.0d0*yhi2+80.0d0*yhi*ylo-40.0d0*ylo2)*sqrryh2)*PINR2
     &    )*f2lo+120.0d0*
     &    (
     &    ((3.0d0*yhi-2.0d0*ylo)*sqrryl2*ylo-sqrryh2*yhi2-2.0d0
     &    *(sqrryh2-sqrryl2)*PINR2)*flo-((2.0d0*yhi-3.0d0*ylo)*sqrryh2*yhi+
     &    sqrryl2*ylo2-2.0d0*(sqrryh2-sqrryl2)*PINR2)*fhi
     &    )+
     &    (
     &    2.0d0*((10.0d0*
     &    yhi2-20.0d0*yhi*ylo+7.0d0*ylo2)*sqrryl2*ylo2+(8.0d0*yhi2-
     &    25.0d0*yhi*ylo+20.0d0*ylo2)*sqrryh2*yhi2+8.0d0*(sqrryh2-sqrryl2)
     &    *r4)+((40.0d0*yhi2-80.0d0*yhi*ylo-43.0d0*ylo2)*sqrryl2-(32.0d0*
     &    yhi2-35.0d0*yhi*ylo-80.0d0*ylo2)*sqrryh2)*PINR2
     &    )*f2hi+15.0d0*
     &    (
     &    (4.0d0*(2.0d0*yhi-ylo)*ylo+3.0d0*PINR2)*
     &    f2lo*yhi-24.0d0*(fhi*ylo-flo*yhi)+
     &    (4.0d0*(yhi-2.0d0*ylo)*yhi-3.0d0*PINR2)*f2hi*ylo
     &    )*asinlo*PINR2

        ans1=-(ans2+ans3)
        ans=ans1/(720.0d0*(yhi-ylo))

        RESULT=RESULT+ans

      ENDDO !IY

      RESULT=2.0D0*RESULT !FAKTOR 2, DA NUR R IN REDUCE ALS BREITE

C NOW CHECK FOR PROBLEMS WITH SPLINES
      SUMP=RESULT

      RETURN
      END
