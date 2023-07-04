*CMZ :  3.08/01 04/04/2019  11.39.33  by  Michael Scheer
*CMZ :  3.07/00 15/03/2019  14.12.22  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.69/00 24/10/2012  16.32.37  by  Michael Scheer
*CMZ :  2.57/03 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.57/00 22/11/2005  10.31.34  by  Michael Scheer
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
      SUBROUTINE circpin_omp(NZ,NY,MZ,MY,SUM,SUMP,ISOUR,kfreq,IWSOUR)
c+seq,gplhint.

c+SELF,IF=F90.
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.
c+SELF.

      use circpinmod

C---  SUBROUTINE TO INTEGRATE SPEC OVER CIRCULAR PINHOLE
C     NZ,NY ARE ARRAY DIMENSIONS
C     INTEGRATION IS DONE FOR MZ AND MY POINTS

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEND.

      INTEGER NZ,NY,IY,IZ,MZ,MY,KLM,KLP,KHM,KHP,IIZ,IIY,ISOUR,kfreq
      INTEGER IWSOUR,IWFREQ,IWSOUR1,IWSOUR2,IWSOUR3,IDUM

      DOUBLE PRECISION SUM,SUMP,DZ,FL,FH,SH,SL,A,B,DY
     &  ,SZ(NDOBSVZP),SUMY(NDOBSVyp),SUMYP(NDOBSVyP)
     &  ,FZ(NDOBSVZP),OBSVZF(NDOBSVZP),ZH,ZL
     &  ,DSUM,DSUML,DSUMH

      save iwsour1,iwfreq,iwsour2

      IF (IF1DIM.NE.0) THEN
        CALL circpin1d_omp(NY,MY,fphir_th,SUM,SUMP)
        RETURN
      ENDIF

      IWSOUR=0

      IF (MZ.LT.2) THEN
        WRITE(6,*) '*** ERROR IN circpin_omp: MZ TOO LOW ***'
        STOP
      ENDIF
      IF (MY.LT.2) THEN
        WRITE(6,*) '*** ERROR IN circpin_omp: MY TOO LOW ***'
        STOP
      ENDIF

C--- LOOP OF VERTICAL DIRECTION

      IIY=0
      DO IY=(NY-MY)/2+1,(NY-MY)/2+MY
        IIY=IIY+1

        SUMY(IIY)=0.0D0
        SUMYP(IIY)=0.0D0

        IF (DABS(OBSVY(IY)-PINCEN(2)).LT.PINR-OBSVDY/4.D0) THEN

C--  SPLINE COEFFICENTS

          IIZ=0
          DO IZ=(NZ-MZ)/2+1,(NZ-MZ)/2+MZ
            IIZ=IIZ+1
            FZ(IIZ)=fphir_th((IY-1)*NZ+IZ)
            OBSVZF(IIZ)=OBSVZ(IZ)
          ENDDO !IZ

          CALL FSPLINDX_omp(OBSVDZ,FZ,IIZ,0.D0,0.D0,SZ)

C--  INTEGRATION LIMITS

          ZH=DSQRT(-(OBSVY(IY)-PINCEN(2))**2+PINR**2)
          ZL=PINCEN(3)-ZH
          ZH=PINCEN(3)+ZH

          CALL SPLINZY_omp(IIZ,ZL,FL,OBSVZF,FZ,SZ,KLM)
          CALL SPLINZY_omp(IIZ,ZH,FH,OBSVZF,FZ,SZ,KHM)

          KLP=KLM+1
          KHP=KHM+1

C--  CONTRIBUTION FORM PARTIAL INTERVALLS

          DZ=OBSVDZ
          DY=OBSVDY

          A=0.0D0
          B=1.0D0

          SH=
     &      -FZ(KLM)*A**2/2.0d0*DZ+FZ(KLP)*B**2/2.0d0*DZ
     &      +DZ**2/6.0d0
     &      *((-DZ*A**4/4.0d0+DZ*A**2/2.)*SZ(KLM)
     &      +(DZ*B**4/4.0d0-DZ*B**2/2.)*SZ(KLP))

          A=(OBSVZF(KLP)-ZL)/OBSVDZ
          B=(ZL-OBSVZF(KLM))/OBSVDZ

          SL=
     &      -FZ(KLM)*A**2/2.0d0*DZ+FZ(KLP)*B**2/2.0d0*DZ
     &      +DZ**2/6.
     &      *((-DZ*A**4/4.0d0+DZ*A**2/2.)*SZ(KLM)
     &      +(DZ*B**4/4.0d0-DZ*B**2/2.)*SZ(KLP))

          DSUML=SH-SL

          IF(
     &        (IWSOUR1.NE.ISOUR.OR.IWFREQ.NE.kfreq)
     &        .AND. ISOUR.GT.0 .AND.
     &        DSUML.LT.0.0) THEN
            IF ( DABS((SH-SL)/(SH+SL)) .GT. 5.D-15 ) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** WARNING IN circpin_omp ***'
              WRITE(LUNGFO,*)
     &          'SPLINE INTEGRATION FAILED FOR LEFT EDGE OF PINHOLE, RESULTS NOT RELIABLE'
              WRITE(LUNGFO,*)'SOURCE POINT AND PHOTON ENERGY:'
     &          ,ISOUR,SNGL(FREQ(kfreq))
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(6,*)
              WRITE(6,*)
              WRITE(6,*)'*** WARNING IN circpin_omp ***'
              WRITE(6,*)
     &          'SPLINE INTEGRATION FAILED FOR LEFT EDGE OF PINHOLE, RESULTS NOT RELIABLE'
              WRITE(6,*)'SOURCE POINT AND PHOTON ENERGY:'
     &          ,ISOUR,SNGL(FREQ(kfreq))
              WRITE(6,*)
              WRITE(6,*)
              IWSOUR1=ISOUR
              IWSOUR=1
              IWFREQ=kfreq
              IW_CIRC=1
            ENDIF    !SH-SL
          ENDIF !IWSOUR1

          A=1.D0
          B=0.D0

          SL=
     &      -FZ(KHM)*A**2/2.0d0*DZ+FZ(KHP)*B**2/2.0d0*DZ
     &      +DZ**2/6.
     &      *((-DZ*A**4/4.0d0+DZ*A**2/2.)*SZ(KHM)
     &      +(DZ*B**4/4.0d0-DZ*B**2/2.)*SZ(KHP))

          A=(OBSVZF(KHP)-ZH)/OBSVDZ
          B=(ZH-OBSVZF(KHM))/OBSVDZ

          SH=
     &      -FZ(KHM)*A**2/2.0d0*DZ+FZ(KHP)*B**2/2.0d0*DZ
     &      +DZ**2/6.
     &      *((-DZ*A**4/4.0d0+DZ*A**2/2.)*SZ(KHM)
     &      +(DZ*B**4/4.0d0-DZ*B**2/2.)*SZ(KHP))


          DSUMH=SH-SL

          IF(
     &        (IWSOUR2.NE.ISOUR.OR.IWFREQ.NE.kfreq)
     &        .AND. ISOUR.GT.0 .AND.
     &        DSUMH.LT.0.0) THEN
            IF ( DABS((SH-SL)/(SH+SL)) .GT. 5.D-15 ) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** WARNING IN circpin_omp ***'
              WRITE(LUNGFO,*)
     &          'SPLINE INTEGRATION FAILED FOR RIGHT EDGE OF PINHOLE, RESULTS NOT RELIABLE'
              WRITE(LUNGFO,*)'SOURCE POINT AND PHOTON ENERGY:'
     &          ,ISOUR,SNGL(FREQ(kfreq))
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(6,*)
              WRITE(6,*)
              WRITE(6,*)'*** WARNING IN circpin_omp ***'
              WRITE(6,*)
     &          'SPLINE INTEGRATION FAILED FOR RIGHT EDGE OF PINHOLE, RESULTS NOT RELIABLE'
              WRITE(6,*)'SOURCE POINT AND PHOTON ENERGY:'
     &          ,ISOUR,SNGL(FREQ(kfreq))
              WRITE(6,*)
              WRITE(6,*)
              IWSOUR2=ISOUR
              IWSOUR=1
              IWFREQ=kfreq
              IW_CIRC=1
            ENDIF
          ENDIF !IWSOUR2

C---  NOW ALL FULL INTERVALLS

          DO IZ=KLP,KHM-1

            DSUM=
     &        OBSVDZ*0.5D0
     &        *(FZ(IZ)+FZ(IZ+1))
     &        -OBSVDZ**3/24.D0
     &        *(SZ(IZ)+SZ(IZ+1))

            IF(
     &          (IWSOUR3.NE.ISOUR.OR.IWFREQ.NE.kfreq)
     &          .AND. ISOUR.GT.0 .AND.
     &          DSUM.LT.0.0) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** WARNING IN circpin_omp ***'
              WRITE(LUNGFO,*)
     &          'SPLINE INTEGRATION FAILED INSIDE PINHOLE, RESULTS NOT RELIABLE'
              WRITE(LUNGFO,*)'SOURCE POINT AND PHOTON ENERGY:'
     &          ,ISOUR,SNGL(FREQ(kfreq))
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(6,*)
              WRITE(6,*)
              WRITE(6,*)'*** WARNING IN circpin_omp ***'
              WRITE(6,*)
     &          'SPLINE INTEGRATION FAILED INSIDE, RESULTS NOT RELIABLE'
              WRITE(6,*)'SOURCE POINT AND PHOTON ENERGY:'
     &          ,ISOUR,SNGL(FREQ(kfreq))
              WRITE(6,*)
              WRITE(6,*)
              IWSOUR3=ISOUR
              IWSOUR=1
              IWFREQ=kfreq
              IW_CIRC=1
              DO IDUM=1,IIY
                SUMYP(IIY)=SUMY(IIY)
              ENDDO
            ENDIF !IWSOUR3

            SUMY(IIY)=SUMY(IIY)+DSUM
            IF (IW_CIRC.NE.0) SUMYP(IIY)=SUMYP(IIY)+ABS(DSUM)

          ENDDO !IZ

          SUMY(IIY)=(SUMY(IIY)+DSUML+DSUMH)

          IF (ZH.GT.0.0D0) THEN
            SUMY(IIY)=SUMY(IIY)/(2.0D0*ZH)
          ELSE
            IF (SUMY(IIY).NE.0.0D0) THEN
              PRINT*,
     &          '*** WARNING IN PINCIRC: INTEGRATION WIDTH=0 BUT CONTRIBUTION NOT'
              PRINT*,'CONTRIBUTION:',SUMYP(IIY)
            ENDIF
          ENDIF

          IF (IW_CIRC.NE.0) THEN
            SUMYP(IIY)=SUMYP(IIY)+ABS(DSUML)+ABS(DSUMH)
            IF (ZH.GT.0.0D0) THEN
              SUMYP(IIY)=SUMYP(IIY)/(2.0D0*ZH)
            ELSE
              IF (SUMYP(IIY).NE.0.0D0) THEN
                PRINT*,
     &            '*** WARNING IN PINCIRC: INTEGRATION WIDTH=0 BUT CONTRIBUTION NOT'
                PRINT*,'CONTRIBUTION:',SUMYP(IIY)
              ENDIF
            ENDIF
          ENDIF

        ELSE !INSIDE PINHOLE

          IZ=NZ/2+1
          SUMY(IIY)=fphir_th((IY-1)*NZ+IZ) ! FOR circpin1d_omp

        ENDIF !INSIDE PINHOLE

      ENDDO !IY

C--- INTEGRATION OVER VERTICAL DIRECTION

      CALL circpin1D_omp(NY,MY,SUMY,SUM,SUMP)

      IF (IW_CIRC.NE.0) THEN
        CALL circpin1d_omp(NY,MY,SUMYP,SUMP,SUMP)
      ENDIF

      RETURN
      END
