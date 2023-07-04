*CMZ :  2.66/07 10/12/2009  12.32.53  by  Michael Scheer
*CMZ :  2.57/05 05/12/2006  11.02.07  by  Michael Scheer
*CMZ :  2.33/07 04/05/2001  15.03.35  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.34.14  by  Michael Scheer
*CMZ :  1.00/00 04/08/97  12.38.30  by  Michael Scheer
*-- Author :  Michael Scheer
      PROGRAM wave_fit_main

C---  MAIN PROGRAM TO USE MINUIT IN DATA-DRIVEN MODE

      IMPLICIT NONE

      EXTERNAL FCN
      EXTERNAL FUTIL
      DOUBLE PRECISION FUTIL

      INTEGER LREAD,LWRITE,LSAVE

C     UNIT 80 F"UR INPUTFILE RESERVIERT

      DATA LREAD/5/
      DATA LWRITE/6/
      DATA LSAVE/7/

      CALL MINTIO(LREAD,LWRITE,LSAVE)

      OPEN(UNIT=LSAVE,FILE='wave-fit.out',STATUS='UNKNOWN')
      OPEN(UNIT=70,FILE='wave-fit.result',STATUS='NEW',RECL=256)

      CLOSE(70)
      OPEN(UNIT=70,FILE='wave-fit.term',STATUS='NEW')
      CLOSE(70)

      CALL MINUIT(FCN,FUTIL)

      CLOSE(LSAVE)

      STOP
      END

      INCLUDE 'fcn.f'
      INCLUDE 'futil.f'
