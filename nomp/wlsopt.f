*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  12.43.40  by  Michael Scheer
*CMZ :  2.56/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.02.52  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.58.52  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.31.54  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.57.14  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.42  by  Michael Scheer
*-- Author : Michael Scheer
C*****************************************************************
      SUBROUTINE WLSOPT(IWLS,IFLAG,E,B0MIN,B0MAX,DB0,EMIKRIT,TAUKRIT,POLKRIT,
     &                    XLAM0MN,XLAM0MX,DXLAM0,FASYMMN,FASYMMX,DFASYM,
     &                    BETA0,RHODIP,DBHOM,TAU0E5,DI2RING,DI5RING,
     &                    ZMAXKRIT,ZMINKRIT,IMODE,LUN,IWLSADI,DISP0,DXHOM)
C****************************************************************
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

C   S/R WLSOPT SUCHT DURCH VARIATION DES B-FELDES, DER LAENGE UND DER
C   ASYMMETRIE DEN OPTIMALEN WLS VOM ASYMM. HALBACH-TYP
C
C   BEDEUTUNG EINIGER VARIABLER:
C
C     IWLS  .NE. 0  PROGRAMM WIRD VON WAVE.FOR GERUFEN
C     IFLAG .NE. 0  MINIMALE EMITTANZ BEZIEHT SICH AUF Emin
C                   ZUR ZEIT NUR MOEGLICH, FALLS DISP0=0
C     IMODE  =   1  SUCHE WLS FUER MINIMALE ENERGIE
C            =   2  SUCHE BESTEN AEQUIVALENTEN AHW

      IMPLICIT NONE

      DOUBLE PRECISION EMI,EMIOPT,EMIRING,EMIRINGK,EMITOTEPK,EMIOPTEP,EMITOPTEP,
     &         EMIRINGEP,EMIRINGEPK,EMITOT,EMITOPT,EMIKRIT
      DOUBLE PRECISION EMIK,EMIOPTK,EMITOTK,EMITOPTK,EMITOTEP,EMIEP,EMIEPK
      DOUBLE PRECISION B0,B0MAX,B0MIN,DB0,B0K,DBHOM,XHOMK,XHOM
      DOUBLE PRECISION RHODIP,RHO0K,GMOM
      DOUBLE PRECISION FASYM,FASYMMN,FASYMMX,DFASYM,FASYMK
      DOUBLE PRECISION E,F,CD,C2,C5,E0,EMINEM,EMINPOL,EMINEP,EMIN,EMINK,GAMMA
      DOUBLE PRECISION ZMAX,ZMAXK,ZMAXEP,ZMAXEPK
      DOUBLE PRECISION BETA0,BETOPT,BETOPTK,F0,F0K,XLP,FB0M,FBETP,B0P,F0P
      DOUBLE PRECISION XLAM0,XLAM0MN,XLAM0MX,DXLAM0,XLAM0K,XKK,XLK,XL2K,DXHOM
      DOUBLE PRECISION TAU1GEV,TAU0E5,TAUMAX,TAU1GEVK,TAUKRIT,TAU,TAUK,TAUEP,TAUEPK
      DOUBLE PRECISION POLLEV,POLLV1G,POLLV1GK,POLLEVK,POLKRIT,POLLEVEP,POLLEVEPK
      DOUBLE PRECISION DI2RING,DI5RING,B3AWLS,DI2WLS,DI2WLSEP,DI2WLSEPK,DI5WLSK
      DOUBLE PRECISION DI5WLSOPT,DI5WLSOPTK,CHI2,CHI2MIN,BETIN,DI5IN,DI2IN  !060891
      DOUBLE PRECISION DI5WLS,DI5WLSEP,DI2ADI,DI5ADI,CHI2ADI,BETOADI,DI5WLSOPTEP
      DOUBLE PRECISION TAUIN,POLLEVIN,ZMAXKRIT,ZMINKRIT,EMIIN,EMIOPTIN
     &                 ,DISP0,DISPOPT,EMI2OPT,EMIT2OPT,BET2OPT,BETUNI,
     &                  BET2UNI,DI2WLSK,
     &                  DISPOPTEP,EMI2OPTEP,EMIT2OPTEP,BET2OPTEP
     &                 ,BETUNIEP,BET2UNIEP,DI5WLS2OPTEP,BETOPTEP,
     &                  F2,F20,EMIT2OPTK,EMI2OPTK,DI5WLS2OPT,BET2OPTK,
     &                  BETUNIK,BET2UNIK,F20K,DI5WLS2OPTK,DISPOPTK,
     &                  BETOPTEPK,DISPOPTEPK,BET2OPTEPK,BETUNIEPK,
     &                  BET2UNIEPK

      INTEGER IB0,NB0,IFA,NFA,IXL,NXL,LUN,IFLAG,IEMIN0,IFOUND,IWLS,IMODE,
     &          IWLSADI

*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'PROGRAMM WLSOPT'
      WRITE(LUN,*)'==============='
      WRITE(LUN,*)


C--- WENN SICH DIE STRAHLENERGIE AENDERT, AENDERT SICH DIE SKALIERTE
C    EXTERNE DISPERSION. DESHALB MUESSTE FUER DISP0.NE.0 UND IFLAG.NE.0
C    DAS PROGRAMM ERWEITERT WERDEN (!? 21.7.92)

      IF(IFLAG.NE.0.AND.DISP0.NE.0.0) THEN
         WRITE(LUN,*)
         WRITE(LUN,*)'*** DISP0.NE.0 UND IEMICRIT.NE.0 ***'
         WRITE(LUN,*)
     &'FUER DISP0.NE.0 MUSS SICH DIE MAX. ZULAESSIGE AENDERUNG AUF E UND NICHT AUF Emin BEZIEHEN!'
         STOP
      ENDIF

      IF(IMODE.EQ.2) THEN
        IF(IWLS.EQ.0) THEN
          WRITE(6,*)  !060891
CC       WRITE(6,*)'BETIN,DI5IN,DI2IN:'   !060891
          WRITE(6,*)'BETIN,DI5IN,DI2IN,TAUIN,POLLEVIN:'   !060891
          READ(5,*)BETIN,DI5IN,DI2IN                !060891
     &      ,TAUIN,POLLEVIN      !27.1.92
        ELSE
          READ(39,*)BETIN,DI5IN,DI2IN        !070891
     &      ,TAUIN,POLLEVIN,EMIIN,EMIOPTIN  !27.1.92
        ENDIF
      ENDIF !IMODE.EQ.2

      GAMMA=E/EMASSG1
      EMIN=E
      GMOM=DSQRT((E-EMASSG1)*(E+EMASSG1))

      IF (LUN.NE.6.AND.LUN.NE.16)
     &  OPEN(UNIT=LUN,FILE='WAVE_WLSOPT.DAT',STATUS='NEW',FORM='FORMATTED')

      CD=EMIKRIT*DI5RING/DI2RING

      NB0=NINT((B0MAX-B0MIN)    /DB0   )
      NXL=NINT((XLAM0MX-XLAM0MN)/DXLAM0)
      NFA=NINT((FASYMMX-FASYMMN)/DFASYM)

C- VARIIERE DIE PARAMETER DES AHW, UM IDEALES GERAET ZU FINDEN

      DO IFA=0,NFA
        FASYM=FASYMMN+IFA*DFASYM
        DO IB0=0,NB0
          B0=B0MIN+IB0*DB0
          DO IXL=0,NXL

            XLAM0=XLAM0MN+IXL*DXLAM0
            XHOM=DSQRT(2.*DBHOM)/(2.* PI1/XLAM0)

C- BERECHNE EMITTANZ UND POLARISATIONSZEIT DES AHW

            CALL EMIT(B0,XLAM0,FASYM,
     &        E,RHODIP,TAU0E5,BETA0,DI2RING,DI5RING,
     &        EMIRING,EMI,EMIOPT,EMITOT,EMITOPT,B3AWLS,
     &        F,F0,BETOPT,TAU,TAU1GEV,POLLEV,POLLV1G,
     &        ZMAX,DI2WLS,DI5WLS,DI5WLSOPT,0,
     &        DISP0,DISPOPT,EMI2OPT,EMIT2OPT,BET2OPT,BETUNI,
     &        BET2UNI,DI5WLS2OPT)

C- START-WERT FUER NUMERISCHE LOESUNG DER NACHFOLGENDEN GLEICHUNG
C  EMINEM IST ENERGIE,FUER DIE GESAMTEMITTANZ GLEICH DER ZULAESSIGEN
C  GESAMTEMITTANZ IST

            IF(IFLAG.NE.0) THEN !1.4.92

              EMINEM=(F*0.5*(BETA0/BETOPT+BETOPT/BETA0)*
     &          (XLAM0*0.3*B0)**3/CD)**(1./3.)

C010891         C5=F*XLAM0**4/4.*(0.3*B0)**5*BETA0/(2.*BETOPT)*(1.+1./FASYM)
              C5=F*XLAM0**4/4.*(0.3*B0)**5*
     &          0.5*(BETA0/BETOPT+BETOPT/BETA0)*(1.+1./FASYM)
              C2=XLAM0/4.*(0.3*B0)**2*(1.+1./FASYM)

C- BESTIMME EMINEM NUMERISCH NACH NEWTON-VERFAHREN. GL. S. LOGBUCH S.194

              CALL EMIN0(CD,C5,C2,DI2RING,DI5RING,E,EMINEM,E0,IEMIN0)

              IF(IEMIN0.NE.1) GOTO 900
              EMINEM=E0

            ENDIF !IFLAG

C- BESTIMME EMINPOL NUMERISCH NACH NEWTON-VERFAHREN. GL. S. LOGBUCH S.193
C  EMINPOL IST DIE ENERGIE, FUER DIE DIE KRITISCHE POLARISATIONSZEIT ERREICHT
C  WIRD

            CALL EMINP(RHODIP,B3AWLS,TAU0E5,TAUKRIT,EMINPOL)

            IF(IFLAG.NE.0) THEN
              EMINEP=DMAX1 (EMINEM,EMINPOL)  !EMITTANZCUT FUER Emin
            ELSE
              EMINEP=EMINPOL !DANN GREIFT EMITTANZCUT FUER E (S.U.)
            ENDIF

C- BERECHNE EMITTANZ UND POL.-ZEIT FUER DIE MINIMALE ENERGIE

            CALL EMIT(B0,XLAM0,FASYM,
     &        EMINEP,RHODIP,TAU0E5,BETA0,DI2RING,DI5RING,
     &        EMIRINGEP,EMIEP,EMIOPTEP,EMITOTEP,EMITOPTEP,
     &        B3AWLS,
     &        F2,F20,BETOPTEP,TAUEP,TAU1GEV,POLLEVEP,POLLV1G,
     &        ZMAXEP,DI2WLSEP,DI5WLSEP,DI5WLSOPTEP,0,
     &        DISP0,DISPOPTEP,EMI2OPTEP,EMIT2OPTEP,BET2OPTEP
     &        ,BETUNIEP,BET2UNIEP,DI5WLS2OPTEP)

            TAUMAX=TAUEP

C--- KANDIDAT GEFUNDEN ?

            IF(IMODE.EQ.2) !260291
     &        CHI2=0.0
     &        +((DI2WLS   -DI2IN)/DI2IN)**2
     &        +((EMI-EMIIN)      /EMIIN)**2
C     &          +((BETOPT   -BETIN)/BETIN)**2  !060891
C     &          +((DI5WLSOPT-DI5IN)/DI5IN)**2
C     &      CHI2=((BETOPT   -BETIN)/BETIN)**2  !060891
C     &          +((DI5WLSOPT/DI2WLS-DI5IN/DI2IN)/(DI5IN/DI2IN))**2
C     &          +((DI2WLS   -DI2IN)/DI2IN)**2
C     &     +(((POLLEV-POLLEVIN)/POLLEVIN)**2     !27.1.92
C     &          +((TAU-TAUIN)/TAUIN)**2)/2.           !27.1.92

            IF(IMODE.EQ.2.AND.CHI2.LT.CHI2MIN) THEN
              IFOUND=1
              EMIN=EMINEP
              EMINK=EMINEP
              EMIRINGK=EMIRING
              EMIRINGEPK=EMIRINGEP
              EMIK=EMI
              EMIEPK=EMIEP
              EMIOPTK=EMIOPT
              EMITOTK=EMITOT
              EMITOTEPK=EMITOTEP
              EMITOPTK=EMITOPT
              DI5WLSK=DI5WLS
              DI5WLSOPTK=DI5WLSOPT !060891

              EMI2OPTK=EMI2OPT
              EMIT2OPTK=EMIT2OPT
              DI5WLS2OPTK=DI5WLS2OPT
              DISPOPTK=DISPOPT
              DISPOPTEPK=DISPOPTEP
              BETOPTEPK=BETOPTEP
              BET2OPTEPK=BET2OPTEP
              BET2OPTK=BET2OPT
              BETUNIEPK=BETUNIEP
              BET2UNIEPK=BET2UNIEP
              BETUNIK=BETUNI
              BET2UNIK=BET2UNI
              F20K=F20
              ZMAXEPK=ZMAXEP
              DI2WLSEPK=DI2WLSEP

              DI2WLSK=DI2WLS
              B0K=B0
              RHO0K=1./(CLIGHT1*B0K/GMOM)
              FASYMK=FASYM
              XLAM0K=XLAM0
              XLK=XLAM0K/2.*(1.+FASYM)
              XL2K=XLK/2.
              XKK=2.* PI1/XLAM0K
              XHOMK=XHOM
              TAUK=TAU
              TAU1GEVK=TAU1GEV
              TAUEPK=TAUEP
              POLLEVEPK=POLLEVEP
              POLLEVK=POLLEV
              POLLV1GK=POLLV1G
              BETOPTK=BETOPT
              F0K=F0
              ZMAXK=ZMAX
              CHI2MIN=CHI2   !060891
            ENDIF !IMODE.EQ.2

            IF(IMODE.EQ.1 .AND.
     &          ZMAX.LE. ZMAXKRIT.and.
     &          ZMAX.GE. ZMINKRIT.and.
     &          XHOM.GE.DXHOM.AND.
     &          TAUMAX.LE.TAUKRIT*1.000001
     &          .AND.
     &          POLLEVEP.GE.POLKRIT
     &          .AND.
     &          (IFLAG.EQ.0.AND.EMITOT.LE.EMIKRIT*EMIRING
     &          .OR.IFLAG.NE.0.AND.EMINEM.LE.E) !DANN IST EMITTANZCUT FU
     &          .AND.
     &          EMINEP.LT.EMIN) THEN
              IFOUND=1
              EMIN=EMINEP
              EMINK=EMINEP
              EMIRINGK=EMIRING
              EMIRINGEPK=EMIRINGEP
              EMIK=EMI
              EMIEPK=EMIEP
              EMIOPTK=EMIOPT
              EMITOTK=EMITOT
              EMITOTEPK=EMITOTEP
              EMITOPTK=EMITOPT
              DI5WLSK=DI5WLS
              DI5WLSOPTK=DI5WLSOPT !060891
              DI2WLSEPK=DI2WLSEP
              DI2WLSK=DI2WLS

              EMI2OPTK=EMI2OPT
              EMIT2OPTK=EMIT2OPT
              DI5WLS2OPTK=DI5WLS2OPT
              DISPOPTK=DISPOPT
              DISPOPTEPK=DISPOPTEP
              BETOPTEPK=BETOPTEP
              BET2OPTK=BET2OPT
              BET2OPTEPK=BET2OPTEP
              BETUNIK=BETUNI
              BETUNIEPK=BETUNIEP
              BET2UNIEPK=BET2UNIEP
              BET2UNIK=BET2UNI
              F20K=F20
              ZMAXEPK=ZMAXEP

              B0K=B0
              RHO0K=1./(CLIGHT1*1.D-9*B0K/GMOM)
              FASYMK=FASYM
              XLAM0K=XLAM0
              XLK=XLAM0K/2.*(1.+FASYM)
              XL2K=XLK/2.
              XKK=2.* PI1/XLAM0K
              XHOMK=XHOM
              TAUK=TAU
              TAU1GEVK=TAU1GEV
              TAUEPK=TAUEP
              POLLEVEPK=POLLEVEP
              POLLEVK=POLLEV
              POLLV1GK=POLLV1G
              BETOPTK=BETOPT
              F0K=F0
              ZMAXK=ZMAX
              CHI2MIN=CHI2   !060891
            ENDIF

900         continue
          ENDDO !FASYM
        ENDDO !B0
      ENDDO !LAMBDA0

      IF (IFOUND.NE.1) THEN
        WRITE(LUN,*) '*** KEINEN PASSENDEN WLS GEFUNDEN ***'
        WRITE(6,*) '*** KEINEN PASSENDEN WLS GEFUNDEN ***'
        STOP
      ENDIF

      IF(IMODE.EQ.2) THEN

        IF(DABS(B0K-B0MIN)     .LT.1.D-3 .AND. B0MIN.NE.B0MAX) THEN
          WRITE(LUN,*)
          WRITE(LUN,*) '*** S/R WLSOPT: B0MIN-LIMIT REACHED ***'
          WRITE(LUN,*)
        ENDIF
        IF(DABS(B0K-B0MAX)     .LT.1.D-3 .AND. B0MIN.NE.B0MAX) THEN
          WRITE(LUN,*)
          WRITE(LUN,*)'*** S/R WLSOPT: B0MAX-LIMIT REACHED ***'
          WRITE(LUN,*)
        ENDIF
        IF(DABS(FASYMK-FASYMMN).LT.1.D-3 .AND. FASYMMN.NE.FASYMMX) THEN
          WRITE(LUN,*)
          WRITE(LUN,*)'*** S/R WLSOPT: FASYMMN-LIMIT REACHED ***'
          WRITE(6,*)'*** S/R WLSOPT: FASYMMN-LIMIT REACHED ***'
          WRITE(LUN,*)
        ENDIF
        IF(DABS(FASYMK-FASYMMX).LT.1.D-3 .AND. FASYMMN.NE.FASYMMX)  THEN
          WRITE(LUN,*)
          WRITE(LUN,*)'*** S/R WLSOPT: FASYMMX-LIMIT REACHED ***'
          WRITE(6,*)'*** S/R WLSOPT: FASYMMX-LIMIT REACHED ***'
          WRITE(LUN,*)
        ENDIF
        IF(DABS(XLAM0K-XLAM0MN).LT.1.D-3 .AND. XLAM0MN.NE.XLAM0MX) THEN
          WRITE(LUN,*)
          WRITE(LUN,*)'*** S/R WLSOPT: XLAM0MN-LIMIT REACHED ***'
          WRITE(6,*)'*** S/R WLSOPT: XLAM0MN-LIMIT REACHED ***'
          WRITE(LUN,*)
        ENDIF
        IF(DABS(XLAM0K-XLAM0MX).LT.1.D-3 .AND. XLAM0MN.NE.XLAM0MX) THEN
          WRITE(LUN,*)
          WRITE(LUN,*)'*** S/R WLSOPT: XLAM0MX-LIMIT REACHED ***'
          WRITE(6,*)'*** S/R WLSOPT: XLAM0MX-LIMIT REACHED ***'
          WRITE(LUN,*)
        ENDIF
      ENDIF
C---------------------------------------------------------------------
      IF(IWLSADI.NE.0)
     &  CALL ADI(GAMMA,B0K,XLAM0K,FASYMK,F0K,DI2WLSEPK,DI5WLSOPTK,
     &  BETOPTK,B0P,XLP,BETOADI,FB0M,F0P,FBETP,
     &  DI2ADI,DI5ADI,CHI2ADI)
C---------------------------------------------------------------------

      IF(
     &      DABS(B0K-B0MIN)     .LT.1.D-3.AND.B0MAX.NE.B0MIN.OR.
     &      DABS(B0K-B0MAX)     .LT.1.D-3.AND.B0MAX.NE.B0MIN.OR.
     &      DABS(XLAM0K-XLAM0MN).LT.1.D-3.AND.XLAM0MX.NE.XLAM0MN.OR.
     &      DABS(XLAM0K-XLAM0MX).LT.1.D-3.AND.XLAM0MX.NE.XLAM0MN.OR.
     &      DABS(FASYMK-FASYMMN).LT.1.D-3.AND.FASYMMX.NE.FASYMMN.OR.
     &      DABS(FASYMK-FASYMMX).LT.1.D-3.AND.FASYMMX.NE.FASYMMN)
     &   THEN
         WRITE(LUN,*)
         WRITE(LUN,*) '*** PARAMETER LIMIT REACHED ***'
         WRITE(LUN,*)
      ENDIF
      WRITE(LUN,*)
      WRITE(LUN,*)'VORGEGEBENE PARAMETER:'
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'ENERGIE, GAMMA:                 ',SNGL(E),SNGL(GAMMA)
      WRITE(LUN,*)'VORGEGEBENE BETAFUNKTION:      ',SNGL(BETA0)
      WRITE(LUN,*)'VORGEGEBENE DISPERSION:        ',SNGL(DISP0)
      WRITE(LUN,*)
      WRITE(LUN,*)'B0min, B0max, dB0:              ',
     &               SNGL(B0MIN),SNGL(B0MAX),SNGL(DB0)
      WRITE(LUN,*)'FASYMMN,FASYMMX,DFASYM:         ',
     &               SNGL(FASYMMN),SNGL(FASYMMX),SNGL(DFASYM)
      WRITE(LUN,*)'LAMBD0min, LAMBDA0max, dLAMBDA0:',
     &           SNGL(XLAM0MN),SNGL(XLAM0MX),SNGL(DXLAM0)

            WRITE(LUN,*)
            WRITE(LUN,*)'RADIUS DER DIPOLE, TAU(1GEV) OHNE WLS:',
     &                   SNGL(RHODIP),SNGL(TAU0E5)
            WRITE(LUN,*)'I2(RING), I5(RING):',
     &                   SNGL(DI2RING),SNGL(DI5RING)

      IF(IMODE.EQ.2) THEN

         WRITE(LUN,*)
         WRITE(LUN,*)'BETIN, DI5IN, DI2IN:'  !060891
         WRITE(LUN,*)SNGL(BETIN),SNGL(DI5IN),SNGL(DI2IN) !060891
         WRITE(LUN,*)'TAUIN, POLLEVIN:'  !270192
         WRITE(LUN,*)SNGL(TAUIN),SNGL(POLLEVIN)  !270192

      ELSEIF (IMODE.EQ.1) THEN

            WRITE(LUN,*)
          WRITE(LUN,*)'MAX. EMI.VERSCHLECHTER.(REL.):  ',SNGL(EMIKRIT)
          IF (IFLAG.EQ.0) THEN
         WRITE(LUN,*) '(BEZIEHT SICH AUF E)'
          ELSE
         WRITE(LUN,*) '(BEZIEHT SICH AUF Emin)'
          ENDIF
          WRITE(LUN,*)
          WRITE(LUN,*)'MAX. POLARISATIONSZEIT:         ',SNGL(TAUKRIT)
          WRITE(LUN,*)'MIN. POLARISATIONSGRAD:         ',SNGL(POLKRIT)
          WRITE(LUN,*)
          WRITE(LUN,*)'OBERE GRENZE FUER ABLAGE:       ',SNGL(ZMAXKRIT)
          WRITE(LUN,*)'UNTERE GRENZE FUER ABLAGE:      ',SNGL(ZMINKRIT)

      ELSE

         STOP '*** S/R WLSOPT: IMODE FALSCH ***'

      ENDIF

      WRITE(LUN,*)
      WRITE(LUN,*)'dx (hom):                       ',SNGL(DXHOM)
      WRITE(LUN,*)'dB/B (hom):                     ',SNGL(DBHOM)
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'DATEN DES OPTIMALEN WLS:'
      WRITE(LUN,*)'------------------------'
      WRITE(LUN,*)
      IF(IMODE.EQ.2)
     &   WRITE(LUN,*)'CHI2MIN:            ',SNGL(CHI2MIN)   !060891
      IF(NB0.LT.2) THEN
         WRITE(LUN,*)
         WRITE(LUN,*) '*** WARNUNG: B0 NICHT VARIIERT ***'
      ENDIF
      WRITE(LUN,*)'B0, RHO0:                       ',SNGL(B0K),SNGL(RHO0K)
      IF(NB0.LT.2) THEN
         WRITE(LUN,*)
      ENDIF
      IF(NFA.LT.2) THEN
         WRITE(LUN,*)
         WRITE(LUN,*) '*** WARNUNG: n NICHT VARIIERT ***'
      ENDIF
      WRITE(LUN,*)'ASYMMETRIE:                     ',SNGL(FASYMK)
      IF(NFA.LT.2) THEN
         WRITE(LUN,*)
      ENDIF
      IF(NXL.LT.2) THEN
         WRITE(LUN,*)
         WRITE(LUN,*) '*** WARNUNG: LAMBA0 NICHT VARIIERT ***'
      ENDIF
      WRITE(LUN,*)'LAMBDA0, K:                     ',SNGL(XLAM0K),SNGL(XKK)
      IF(NXL.LT.2) THEN
         WRITE(LUN,*)
      ENDIF
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'GESAMTLAENGE, HALBE LAENGE:     ',SNGL(XLK),SNGL(XL2K)
      WRITE(LUN,*)'2*Xhom, 2*Xhom/RHO0K*GAMMA:     ',
     &              SNGL(2.*XHOMK),SNGL(2.*XHOMK/RHO0K*GAMMA)
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'I2(WLS):                        ',SNGL(DI2WLSK)
      WRITE(LUN,*)'I5(WLS):                        ',SNGL(DI5WLSK)
      WRITE(LUN,*)'I5(WLS) FUER OPT. BETA:         ',SNGL(DI5WLSOPTK)
      WRITE(LUN,*)'I5(WLS) FUER OPT. BETA UND ETA: ',SNGL(DI5WLS2OPTK)
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'EMI(WLS):                       ',SNGL(EMIK)
      WRITE(LUN,*)'EMI(WLS) FUER OPT. BETA:        ',SNGL(EMIOPTK)
      WRITE(LUN,*)'EMI(WLS) FUER OPT. BETA UND ETA:',SNGL(EMI2OPTK)
      WRITE(LUN,*)
      WRITE(LUN,*)'GESAMTEMITTANZ:                       ',SNGL(EMITOTK)
      WRITE(LUN,*)'GESAMTEMITTANZ BEI OPT. BETA:         ',SNGL(EMITOPTK)
      WRITE(LUN,*)'GESAMTEMITTANZ BEI OPT. BETA UND ETA: ',SNGL(EMIT2OPTK)
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'EMITTANZAENDERUNG FUER E:                  ',
     &  SNGL(EMITOTK/EMIRINGK)
      WRITE(LUN,*)'EMITTANZAENDERUNG FUER E UND OPT. BETA:    ',
     &  SNGL(EMITOPTK/EMIRINGK)
      WRITE(LUN,*)'EMITTANZAENDERUNG FUER E,OPT. BETA UND ETA:',
     &  SNGL(EMIT2OPTK/EMIRINGK)
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'OPT. DISPERSION:                       ',SNGL(DISPOPTK)
      WRITE(LUN,*)'OPT. BETAFUNKTION:                     ',SNGL(BETOPTK)
      WRITE(LUN,*)'OPT. BETAFUNKTION FUER OPT. DISPERSION:',SNGL(BET2OPTK)
      WRITE(LUN,*)'NEUTRALE BETAFUNKTION:                 ',SNGL(BETUNIK)
      WRITE(LUN,*)'NEUTRALE BETAFUNKTION FUER OPT. ETA:   ',SNGL(BET2UNIK)
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'ABLAGE:                             ',SNGL(ZMAXK)
      WRITE(LUN,*)'FOKALLAENGE (FALLS 1/Fx=0;NAEH.)    ',SNGL(1./DI2WLSK)
      WRITE(LUN,*)'ABGESTRAHLTE LEISTUNG (kWATT/AMPERE)',
     &               SNGL(14.085*DI2WLSK*E**4)
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'MIN. ENERGIE:                   ',SNGL(EMINK)
      WRITE(LUN,*)'--------------------------------'
      WRITE(LUN,*)
      WRITE(LUN,*)'POLARISATIONSZEIT:              ',SNGL(TAUK)
      WRITE(LUN,*)'POLARISATIONSGRAD:              ',SNGL(POLLEVK)
      WRITE(LUN,*)'POLARISATIONSZEIT (Emin):       ',SNGL(TAUEPK)
      WRITE(LUN,*)'POLARISATIONSGRAD (Emin):       ',SNGL(POLLEVEPK)
      WRITE(LUN,*)'POLARISATIONSZEIT (1 GeV):      ',SNGL(TAU1GEVK)
      WRITE(LUN,*)'POLARISATIONSGRAD (1 GeV):      ',SNGL(POLLV1GK)
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'EMITTANZAENDERUNG FUER Emin:    ',
     &               SNGL(EMITOTEPK/EMIRINGEPK)
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'OPT. DISPERSION FUER Emin:             ',SNGL(DISPOPTEPK)
      WRITE(LUN,*)'OPT. BETAFUNKTION FUER Emin:           ',SNGL(BETOPTEPK)
      WRITE(LUN,*)'OPT. BETAFUNKTION FUER OPT. DISPERSION:',SNGL(BET2OPTEPK)
      WRITE(LUN,*)'NEUTRALE BETAFUNKTION FUER Emin:       ',SNGL(BETUNIEPK)
      WRITE(LUN,*)'NEUTRALE BETAFUNKTION FUER OPT. ETA:   ',SNGL(BET2UNIEPK)
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'ABLAGE FUER Emin:                         ',
     &  SNGL(ZMAXEPK)
      WRITE(LUN,*)'FOKALLAENGE FUER Emin (FALLS 1/Fx=0;NAEH.)',
     &  SNGL(1./DI2WLSEPK)
      WRITE(LUN,*)'ABGESTRAHLTE LEISTUNG (kWATT/AMPERE)      ',
     &               SNGL(14.085*DI2WLSEPK*EMINK**4)


      IF (IWLSADI.NE.0) THEN
      WRITE(LUN,*)'WERTE DES ADI:'
      WRITE(LUN,*)'B0, LPLUS, m:        ',SNGL(B0P),SNGL(XLP),SNGL(FB0M)
      WRITE(LUN,*)'I2(adi),I5opt(adi):  ',SNGL(DI2ADI),SNGL(DI5ADI)
      WRITE(LUN,*)'BETAopt(adi):        ',SNGL(BETOADI)
        WRITE(LUN,*)'CHI2(adi):           ',SNGL(CHI2ADI)
      WRITE(LUN,*)
      ENDIF

      IF (LUN.NE.6.AND.LUN.NE.16)CLOSE(LUN)

      IF (LUN.EQ.16) THEN
C     WRITE(6,*)
C     WRITE(6,*)'PROGRAMM WLSOPT'
C     WRITE(6,*)'==============='
C     WRITE(6,*)
C     WRITE(6,*)
      IF(
     &      DABS(B0K-B0MIN)     .LT.1.D-3.AND.B0MAX.NE.B0MIN.OR.
     &      DABS(B0K-B0MAX)     .LT.1.D-3.AND.B0MAX.NE.B0MIN.OR.
     &      DABS(XLAM0K-XLAM0MN).LT.1.D-3.AND.XLAM0MX.NE.XLAM0MN.OR.
     &      DABS(XLAM0K-XLAM0MX).LT.1.D-3.AND.XLAM0MX.NE.XLAM0MN.OR.
     &      DABS(FASYMK-FASYMMN).LT.1.D-3.AND.FASYMMX.NE.FASYMMN.OR.
     &      DABS(FASYMK-FASYMMX).LT.1.D-3.AND.FASYMMX.NE.FASYMMN)
     &   THEN
         WRITE(6,*)
         WRITE(6,*) '*** PARAMETER LIMIT REACHED ***'
         WRITE(6,*)
      ENDIF
      IF (IFOUND.NE.1) STOP '*** KEINEN PASSENDEN WLS GEFUNDEN ***'
C     WRITE(6,*)
C     WRITE(6,*)'VORGEGEBENE PARAMETER:'
C     WRITE(6,*)
C     WRITE(6,*)'ENERGIE, GAMMA:                 ',SNGL(E),SNGL(GAMMA)
      WRITE(6,*)'B0min, B0max, dB0:              ',
     &               SNGL(B0MIN),SNGL(B0MAX),SNGL(DB0)
      WRITE(6,*)'FASYMMN,FASYMMX,DFASYM:         ',
     &               SNGL(FASYMMN),SNGL(FASYMMX),SNGL(DFASYM)
      WRITE(6,*)'LAMBD0min, LAMBDA0max, dLAMBDA0:',
     &           SNGL(XLAM0MN),SNGL(XLAM0MX),SNGL(DXLAM0)

      IF(IMODE.EQ.1) THEN

         WRITE(6,*)'BETIN, DI5IN, DI2IN:'  !060891
         WRITE(6,*)SNGL(BETIN),SNGL(DI5IN),SNGL(DI2IN) !060891
         WRITE(6,*)

      ELSEIF (IMODE.EQ.2) THEN

          WRITE(6,*)'MAX. EMI.VERSCHLECHTER.(REL.):  ',SNGL(EMIKRIT)
          IF (IFLAG.EQ.0) THEN
         WRITE(6,*) '(BEZIEHT SICH AUF E)'
          ELSE
         WRITE(6,*) '(BEZIEHT SICH AUF Emin)'
          ENDIF
C         WRITE(6,*)'MAX. POLARISATIONSZEIT:         ',SNGL(TAUKRIT)
C         WRITE(6,*)'MIN. POLARISATIONSGRAD:         ',SNGL(POLKRIT)

      ELSE

         STOP '*** S/R WLSOPT: IMODE FALSCH ***'

      ENDIF

C     WRITE(6,*)'dB/B (hom):                     ',SNGL(DBHOM)
C     WRITE(6,*)
C     WRITE(6,*)'DATEN DES OPTIMALEN WLS:'
C     WRITE(6,*)
      IF(IMODE.EQ.2)
     &  WRITE(6,*)'CHI2MIN:            ',SNGL(CHI2MIN)   !060891
      IF(NB0.LT.2) THEN
         WRITE(6,*)
         WRITE(6,*) '*** WARNUNG: B0 NICHT VARIIERT ***'
      ENDIF
      WRITE(6,*)'B0, RHO0:                       ',SNGL(B0K),SNGL(RHO0K)
      IF(NFA.LT.2) THEN
         WRITE(6,*)
         WRITE(6,*) '*** WARNUNG: n NICHT VARIIERT ***'
      ENDIF
      WRITE(6,*)'ASYMMETRIE:                     ',SNGL(FASYMK)
      IF(NXL.LT.2) THEN
         WRITE(6,*)
         WRITE(6,*) '*** WARNUNG: LAMBA0 NICHT VARIIERT ***'
      ENDIF
      WRITE(6,*)'LAMBDA0, K:                     ',SNGL(XLAM0K),SNGL(XKK)
C     WRITE(6,*)'GESAMTLAENGE, HALBE LAENGE:     ',SNGL(XLK),SNGL(XL2K)
      WRITE(6,*)'2*Xhom, 2*Xhom/RHO0K*GAMMA:     ',
     &              SNGL(2.*XHOMK),SNGL(2.*XHOMK/RHO0K*GAMMA)
C     WRITE(6,*)'I2(WLS):                        ',SNGL(DI2WLSEPK)
C     WRITE(6,*)'I5(WLS), OPT. I5(WLS):          ',
C     &               SNGL(DI5WLSK),SNGL(DI5WLSOPTK)
C     WRITE(6,*)'EMI(WLS), OPT. EMI(WLS):        ',SNGL(EMIK),SNGL(EMIOPTK)
C     WRITE(6,*)'EMITTANZ, OPT. EMITTANZ:        ',
C     &               SNGL(EMITOTK),SNGL(EMITOPTK)
C     WRITE(6,*)'EMITTANZAENDERUNG FUER E:       ',SNGL(EMITOTK/EMIRINGK)
C     WRITE(6,*)'EMITTANZAENDERUNG FUER Emin:    ',
C     &               SNGL(EMITOTEPK/EMIRINGEPK)
C     WRITE(6,*)
C     WRITE(6,*)'VORGEGEBENE BETAF., OPT. BETAF.:',
C     &               SNGL(BETA0),SNGL(BETOPTK)
C     WRITE(6,*)'ABLAGE:                         ',SNGL(ZMAXK)
C     WRITE(6,*)'FOKALLAENGE (FALLS 1/Fx=0;NAEH.)',SNGL(1./DI2WLSEPK)
C     WRITE(6,*)
C     WRITE(6,*)'MIN. ENERGIE:                   ',SNGL(EMINK)
C     WRITE(6,*)'POLARISATIONSZEIT:              ',SNGL(TAUK)
C     WRITE(6,*)'POLARISATIONSGRAD:              ',SNGL(POLLEVK)
C     WRITE(6,*)'POLARISATIONSZEIT (Emin):       ',SNGL(TAUEPK)
C     WRITE(6,*)'POLARISATIONSGRAD (Emin):       ',SNGL(POLLEVEPK)
C     WRITE(6,*)'POLARISATIONSZEIT (1 GeV):      ',SNGL(TAU1GEVK)
C     WRITE(6,*)'POLARISATIONSGRAD (1 GeV):      ',SNGL(POLLV1GK)
C     WRITE(6,*)
      IF (IWLSADI.NE.0) THEN
      WRITE(6,*)'WERTE DES ADI:'
      WRITE(6,*)'B0, LPLUS, m:        ',SNGL(B0P),SNGL(XLP),SNGL(FB0M)
      WRITE(6,*)'I2(adi),I5opt(adi):  ',SNGL(DI2ADI),SNGL(DI5ADI)
      WRITE(6,*)'BETAopt(adi):        ',SNGL(BETOADI)
      WRITE(6,*)'CHI2(adi):           ',SNGL(CHI2ADI)
      WRITE(6,*)
      ENDIF
      ENDIF

      RETURN
      END
