*CMZ : 00.00/11 07/06/2011  13.53.13  by  Michael Scheer
*CMZ : 00.00/02 21/06/2004  11.50.19  by  Michael Scheer
*-- Author :    Michael Scheer   16/11/99

      SUBROUTINE UTIL_FLUX_TO_DOSE(E,F,DE,AREA,DOSE,IERR)

C CONVERT FLUX F=dN/dE[eV]/dt[s] IN INTERVALL DE ON AREA[m*m] TO DOSE [mSv/h]
C PHOTON ENERGY E [eV]

      IMPLICIT NONE
      INTEGER N,KLO,KHI,K,IERR,ICAL

      PARAMETER (N=1000)

      DOUBLE PRECISION E,F,DE,Y,X,H,AREA
      DOUBLE PRECISION XA(N),YA(N),BB
      DOUBLE PRECISION DENDUM,ECHARGE1,XA1OLD,XANOLD,YA1OLD,YANOLD,DOSE

      CHARACTER(50) CDUM
      INTEGER IFREQ,KLOLD,NOLD
      INTEGER NMU,IMU

      DATA ICAL/0/
      DATA ECHARGE1/1.602176462D-19/

      IF (ICAL.EQ.0) THEN

        OPEN(UNIT=99,FILE='ABSORPDOSE.RP',STATUS='OLD')

        READ(99,'(A50)') CDUM
        READ(99,*) DENDUM
        READ(99,*) NMU

        IF (NMU.GT.N) STOP '*** UTIL_FLUX_TO_DOSE: DIMENSION EXCEEDED  ***'

        DO IMU=1,NMU
          READ(99,*) XA(IMU),YA(IMU)
        ENDDO !NMU

        ICAL=1

      ENDIF   !ICAL


      X=F*E*ECHARGE1

      IF(XA(1).ge.XA(N)) then
        WRITE(6,*)
        WRITE(6,*) '*** UTIL_FLUX_TO_DOSE: E not increasing ***'
        WRITE(6,*)
        WRITE(6,*) 'E-RANGE:',XA(1),'-',XA(N)
        WRITE(6,*) 'E:      ',X
        WRITE(6,*)
        IERR=1
        stop
      ENDIF

      IF(XA(1).LT.XA(N).AND.(X.LT.XA(1).OR.X.GT.XA(N))
     &    .OR.
     &    XA(N).LT.XA(1).AND.(X.LT.XA(N).OR.X.GT.XA(1))) THEN
        WRITE(6,*)
        WRITE(6,*) '*** UTIL_FLUX_TO_DOSE: E OUT OF RANGE ***'
        WRITE(6,*)
        WRITE(6,*) 'E-RANGE:',XA(1),'-',XA(N)
        WRITE(6,*) 'E:      ',X
        WRITE(6,*)
        IERR=1
        RETURN
      ENDIF

      KLO=1
      KHI=N

1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
        GOTO 1
      ENDIF

      H=XA(KHI)-XA(KLO)

C INTERPOLATION OF Y BY Y=AA*X**BB

      IF (H.NE.0.) THEN
            BB=DLOG(YA(KHI)/YA(KLO))/DLOG(XA(KHI)/XA(KLO))
            Y=YA(KLO)*DEXP(BB*(DLOG(X/XA(KLO))))
      ELSE
C           IF (YA(KLO).NE.YA(KHI)) STOP '*** UTIL_FLUX_TO_DOSE: Bad Input ***'
            Y=YA(KLO)
      ENDIF

      KLOLD=KLO
      NOLD=N
      XA1OLD=XA(1)
      XANOLD=XA(N)
      YA1OLD=YA(1)
      YANOLD=YA(N)

      DOSE=X*Y*DE

      RETURN
      END
