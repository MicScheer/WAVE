*CMZ : 00.00/16 20/08/2014  14.15.38  by  Michael Scheer
*CMZ : 00.00/07 12/05/2010  13.07.27  by  Michael Scheer
*CMZ : 00.00/02 25/03/2002  13.43.17  by  Michael Scheer
*-- Author :    Michael Scheer   04/05/99

      SUBROUTINE UTIL_HARM_K(E,RLAMBD,H,DFLEC,B)

      IMPLICIT NONE

*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

C--- CALCULATES K PARAMETER DFLEC FROM FIRST HARMONIC H[EV] AND PERIODLENGTH
C    RLAMBD [CM]

        REAL*8 E,RLAMBD,H,DFLEC,B,RLAMB1,GAMMA

*KEEP,phycon1.
      include 'phycon1.cmn'
*KEND.

c        DFLEC=2.*(950./H*E*E/RLAMBD-1.)
        RLAMB1=WTOE1/H*1.0d-7
        GAMMA=E/EMASSG1
        DFLEC=2.0d0*(RLAMB1/RLAMBD*2.0d0*gamma**2-1.0d0)
      IF (DFLEC.GT.0.) THEN
          DFLEC=SQRT(DFLEC)
      ELSE
          DFLEC=0.0
      ENDIF

c     B=DFLEC/0.934D0/RLAMBD
        B=DFLEC/(ECHARGE1*RLAMBD/100.0d0/(2.*PI1*EMASSKG1*CLIGHT1))

      RETURN
      END
