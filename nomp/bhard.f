*CMZ :  3.05/11 14/08/2018  12.28.54  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.52/11 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/09 08/03/2000  17.28.37  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.36.07  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.27  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.33  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BHARD(XS,XE,GAMMA,ZI,ZIP,YI,YIP,      ZF,ZFP,YF,YFP)

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

C--- SUBROUTINE BERECHNET ZUNAECHST DIE LINEARE TRANSFERMATRIX EINES HARD-EDGE
C             WLS UND DANN BEI JEDEM AUFRUF DIE AUS DEN KOORDINATEN XI,XIP...
C        DIE ENDKOORDINATEN XF,XFP,...

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,bfeld.
      include 'bfeld.cmn'
*KEND.

      INTEGER ICAL,IMAT,ICOMP,IMAG,JMAG,I,IDRIFT
      DOUBLE PRECISION T3,DL,DL0,RHO,DPHI,XVI,XVF,EMASS,CLIGHT,XI,ZI,XIP,ZIP
     &        ,TWORK,XS,XE,XF,ZF,ZFP,XFP,YF,YFP,YI,YIP,GAMMA,P

        DIMENSION T3(4,4,29) !TRANSFER-MATRIZEN FUER 7 RECHTECKMAG. INCL.
        !WEDGES UND 6 DRIFTS
     &    ,TWORK(4*4*29)

        DIMENSION DL(7),    !LAENGE DER BAHN IM MAGNETEN
     &    DL0(7),   !LAENGE DES MAGNETEN
     &    RHO(7),   !BAHNRADIUS IM MAGNETEN
     &    DPHI(7),  !WINKELAENDERUNG IM MAGNETEN
     &        XVI(4),XVF(4) !HILFSVEKTOREN

      EQUIVALENCE (T3,TWORK)

      DATA ICAL/0/

      DATA EMASS /0.5110034D-3/
      DATA CLIGHT/2.99792458D8/

C--- KOORDINATEN-SYSTEM WECHSELN

      XI=-ZI
      XIP=-ZIP

C---------  ICAL .NE. 1  ----------------------------------------------------
      IF (ICAL.NE.1) THEN  !MATRIX BERECHEN AUS DATEN VON NAMELIST BBFELD

        IF(XS.GT.XM1.OR.XE.LT.XP7)
     &    STOP'*** SR BHARD: XSTART.GT.XM1.OR.XEND.LT.XP7 ***'

C        IN DIESEM FALL NUR 5 MATRIZEN KORREKT
      IF(IKBFORM.NE.0)
     &  STOP '*** STOP SR BHARD: IKBFORM.NE.0 ***'

C--- BERECHNE FUER JEDEN MAGNETEN RHO, LAENGE DER BAHN UND ABLENKWINKEL

          DL0(1)=XP1-XM1
          DL0(2)=XP2-XM2
          DL0(3)=XP3-XM3
          DL0(4)=XP4-XM4
          DL0(5)=XP5-XM5
          DL0(6)=XP6-XM6
          DL0(7)=XP7-XM7

          P=GAMMA*EMASS*1.D9    !IMPULS IN eV


          DO I=1,7
         RHO(I)=0.D0
            ENDDO

          IF(BBY1.NE.0.D0) RHO(1)=P/(CLIGHT*BBY1) !RHO KANN NEGATIV SEIN
          IF(BBY2.NE.0.D0) RHO(2)=P/(CLIGHT*BBY2)
          IF(BBY3.NE.0.D0) RHO(3)=P/(CLIGHT*BBY3)
          IF(BBY4.NE.0.D0) RHO(4)=P/(CLIGHT*BBY4)
          IF(BBY5.NE.0.D0) RHO(5)=P/(CLIGHT*BBY5)
          IF(BBY6.NE.0.D0) RHO(6)=P/(CLIGHT*BBY6)
          IF(BBY7.NE.0.D0) RHO(7)=P/(CLIGHT*BBY7)

          DO I=1,4*4*25
         TWORK(I)=0.D0
          ENDDO

C--- DRIFTS
          DO IDRIFT=1,29,4
         T3(1,1,IDRIFT)=1.D0
         T3(2,2,IDRIFT)=1.D0
         T3(3,3,IDRIFT)=1.D0
         T3(4,4,IDRIFT)=1.D0
          END DO

          T3(1,2,1) = XM1-XS
          T3(3,4,1) = XM1-XS

          T3(1,2,5) = XM2-XP1
          T3(3,4,5) = XM2-XP1
          T3(1,2,9) = XM3-XP2
          T3(3,4,9) = XM3-XP2
          T3(1,2,13)= XM4-XP3
          T3(3,4,13)= XM4-XP3
          T3(1,2,17)= XM5-XP4
          T3(3,4,17)= XM5-XP4
          T3(1,2,21)= XM6-XP5
          T3(3,4,21)= XM6-XP5
          T3(1,2,25)= XM7-XP6
          T3(3,4,25)= XM7-XP6

          T3(1,2,29)= XE-XP7
          T3(3,4,29)= XE-XP7

C--- MAGNETE (SECTOR)

          DO I=1,7
            IF(RHO(I).NE.0.D0) THEN
              DPHI(I)=2.D0*DASIN(DL0(I)/(2.D0*RHO(I)))
              DL(I)=RHO(I)*DPHI(I)
            ELSE
              DL(I)=DL0(I)
              DPHI(I)=0.D0
            ENDIF
          ENDDO

          JMAG=0
          DO IMAG=3,27,4
             JMAG=JMAG+1
             T3(1,1,IMAG)=DCOS(DPHI(JMAG))
             T3(2,2,IMAG)=T3(1,1,IMAG)
             T3(1,2,IMAG)=DL0(JMAG)

             IF(RHO(JMAG).NE.0.0) THEN
               T3(1,2,IMAG)=RHO(JMAG)*DSIN(DPHI(JMAG))
               T3(2,1,IMAG)=-DSIN(DPHI(JMAG))/RHO(JMAG)
             ENDIF

             T3(3,3,IMAG)=1.D0
             T3(4,4,IMAG)=1.D0
             T3(3,4,IMAG)=DL(JMAG)
          ENDDO

C--- WEDGES; DER WINKEL DES KEILES WIRD UEBER DEN EINGELESENEN FAKTOR AUS DPHI
C            BERECHNET

          JMAG=0
          DO I=3,27,4
             JMAG=JMAG+1
             T3(1,1,I-1)=1.D0
             T3(2,2,I-1)=1.D0
             T3(3,3,I-1)=1.D0
             T3(4,4,I-1)=1.D0
             T3(1,1,I+1)=1.D0
             T3(2,2,I+1)=1.D0
             T3(3,3,I+1)=1.D0
             T3(4,4,I+1)=1.D0

          IF (RHO(JMAG).NE.0.0)
     &        T3(2,1,I-1)= DTAN(WEDFACL(JMAG)*DPHI(JMAG))/RHO(JMAG)
          IF (RHO(JMAG).NE.0.0)
     &        T3(2,1,I+1)= DTAN(WEDFACR(JMAG)*DPHI(JMAG))/RHO(JMAG)
          IF (RHO(JMAG).NE.0.0)
     &        T3(4,3,I-1)=-DTAN(WEDFACL(JMAG)*DPHI(JMAG))/RHO(JMAG)
          IF (RHO(JMAG).NE.0.0)
     &        T3(4,3,I+1)=-DTAN(WEDFACR(JMAG)*DPHI(JMAG))/RHO(JMAG)
          ENDDO


      END IF
C---------  ICAL .NE. 1  ----------------------------------------------------

C--- XF,XFP... BERECHNEN

      XVI(1)=XI
      XVI(2)=XIP
      XVI(3)=YI
      XVI(4)=YIP

      DO IMAT=1,29
        DO ICOMP=1,4
          XVF(ICOMP)=
     &      T3(ICOMP,1,IMAT)*XVI(1)+T3(ICOMP,2,IMAT)*XVI(2)+
     &      T3(ICOMP,3,IMAT)*XVI(3)+T3(ICOMP,4,IMAT)*XVI(4)
        END DO
        DO ICOMP=1,4
          XVI(ICOMP)=XVF(ICOMP)
        ENDDO
      END DO

      XF=XVF(1)
      XFP=XVF(2)
      YF=XVF(3)
      YFP=XVF(4)

C--- KOORDINATEN-SYSTEM WECHSELN

      ZF=-XF
      ZFP=-XFP

      ICAL=1

      RETURN
      END
