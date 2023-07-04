*CMZ :  3.05/04 04/07/2018  15.30.22  by  Michael Scheer
*CMZ :  3.02/04 02/12/2014  16.21.37  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.66/20 22/11/2011  10.34.00  by  Michael Scheer
*CMZ :  2.61/02 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  10.41.29  by  Michael Scheer
*CMZ :  2.52/16 21/01/2005  17.08.15  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  10.56.48  by  Michael Scheer
*CMZ :  2.20/01 29/11/2000  18.34.26  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.12/00 27/05/99  18.01.14  by  Michael Scheer
*CMZ :  2.10/01 07/05/99  12.21.34  by  Michael Scheer
*CMZ : 00.01/10 20/08/96  12.00.15  by  Michael Scheer
*CMZ : 00.01/08 21/06/95  17.17.07  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.22.06  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.48.13  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.17  by  Michael Scheer
*-- Author : Michael Scheer
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
C*******************************************************************************
      Subroutine BMOVE(R0X,R0Y,R0Z,V0X,V0Y,V0Z,BXIN,BYIN,BZIN,DELTAT,
     &             X,Y,Z,VX,VY,VZ,VXP,VYP,VZP,Gamma,ICHARGE,BMOVECUT,IUSTEP)

C
C     SEE ALSO MYBMOVE
C
C
C*******************************************************************************
C
C BERECHNET DIE EXAKTE 3-DIM Trajektorie EINES TEILCHENS IN EINEM
C 3-DIM HOMOGENEN Magnetfeld
C  NEUESTE VERSION VOM 15.5.1985
C
C LAENGEN IN METERN
C GeschwindigkeitEN IN M/SEC
C ZEIT IN SEKUNDEN
C B-FELDER IN TESLA V SEC/M**2
C
C*******************************************************************

C20.8.96      ImpLicit DOUBLE PRECISION(A-Z)

        IMPLICIT NONE

        DOUBLE PRECISION BX,BY,BZ,BSQ,BBET,BUX,BUY,BUZ,V0SQ,V0X,V0Y,V0Z,
     &         V0BET,VPAR,VPARX,VPARY,VPARZ,VSENK,DELTAT,BXIN,BYIN,BZIN
        DOUBLE PRECISION N1X,N1Y,N1Z,N2X,N2Y,N2Z,RHO  !DAS WAR ICH NICHT

        DOUBLE PRECISION X,Y,Z,VX,VY,VZ,R0X,R0Y,R0Z,VXP,VYP,VZP,ZYK,GAMMA,SZ,CZ
        DOUBLE PRECISION BMOVECUT,dgammao

        double precision dgamma,acc

      INTEGER ICHARGE,IUSTEP

C      DATA QEL,CLIGHT,EMASS/1.6021892D-19,2.99792458D08, 9.109534D-31/

*KEEP,efield.
      include 'efield.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.
C
C KOMPONENTEN DER EINHEITSVEKTOREN IN
C B-RICHTUNG
C

      IF (ICHARGE.GT.0) THEN
          BX=-BXIN
          BY=-BYIN
          BZ=-BZIN
      ELSE
          BX=BXIN
          BY=BYIN
          BZ=BZIN
      ENDIF

        IF(DABS(BX).LT.BMOVECUT.AND.DABS(BY).LT.BMOVECUT
     &    .AND.DABS(BZ).LT.BMOVECUT) GO TO 30

      BSQ=BX*BX+BY*BY+BZ*BZ
      BBET=SQRT(BSQ)
      BUX=BX/BBET
      BUY=BY/BBET
      BUZ=BZ/BBET

C
C   BETRAG VON V0 PARALELL UND SENKRECHT
C
      V0SQ=V0X*V0X+V0Y*V0Y+V0Z*V0Z
      V0BET=SQRT(V0SQ)
      VPAR=(V0X*BX+V0Y*BY+V0Z*BZ)/BBET
      VSENK=V0BET-VPAR
      IF (VSENK.LT.0.0) THEN
          VSENK=0.0
      ELSE
          VSENK=SQRT((V0BET+VPAR)*VSENK)
      ENDIF

C
C   VEKTOR N1 BERECHNEN
C
      IF(VSENK.EQ.0.) N1X=0.
      IF(VSENK.EQ.0.) N1Y=0.
      IF(VSENK.EQ.0.) N1Z=0.
      IF(VSENK.EQ.0.) GO TO 10

      N1X=(V0X-VPAR*BUX)/VSENK
      N1Y=(V0Y-VPAR*BUY)/VSENK
      N1Z=(V0Z-VPAR*BUZ)/VSENK

10    CONTINUE

C
C  VEKTOR N2=(BUX,BUY,BUZ) KREUZ N1
C
      N2X = BUY*N1Z - BUZ*N1Y
      N2Y = BUZ*N1X - BUX*N1Z
      N2Z = BUX*N1Y - BUY*N1X

C
C  VPAR
C
      VPARX=VPAR*BUX
      VPARY=VPAR*BUY
      VPARZ=VPAR*BUZ

C
C ZYKLOTRONFREQUENZ
C
      ZYK= (ECHARGE1/(Gamma*EMASSKG1))*BBET

C
C
      SZ=DSIN(ZYK*DELTAT)
      CZ=DCOS(ZYK*DELTAT)

C
C  ZEITABLEITUNG DER Geschwindigkeit
C
c This is at the beginning of the step
cerror 4.7.2018      VXP=VSENK*ZYK*N2X
cerror 4.7.2018      VYP=VSENK*ZYK*N2Y
cerror 4.7.2018      VZP=VSENK*ZYK*N2Z

C
C V(DELTAT),X(DELTAT) BERECHNEN
C
      VX=VPARX + VSENK*(N1X*CZ+N2X*SZ)
      VY=VPARY + VSENK*(N1Y*CZ+N2Y*SZ)
      VZ=VPARZ + VSENK*(N1Z*CZ+N2Z*SZ)

      IF(VSENK.EQ.0.) X=R0X+VPARX*DELTAT
      IF(VSENK.EQ.0.) Y=R0Y+VPARY*DELTAT
      IF(VSENK.EQ.0.) Z=R0Z+VPARZ*DELTAT
      IF(VSENK.EQ.0.) GO TO 20

C
C X(DELTAT) BERECHNEN
C

      RHO=VSENK/ZYK

      X=R0X+RHO*(N2X+N1X*SZ-N2X*CZ)+VPARX*DELTAT
      Y=R0Y+RHO*(N2Y+N1Y*SZ-N2Y*CZ)+VPARY*DELTAT
      Z=R0Z+RHO*(N2Z+N1Z*SZ-N2Z*CZ)+VPARZ*DELTAT

20    CONTINUE
      GOTO 999

30    CONTINUE

      VXP=0.D0
      VYP=0.D0
      VZP=0.D0

      X=R0X+V0X*DELTAT
      Y=R0Y+V0Y*DELTAT
      Z=R0Z+V0Z*DELTAT

      VX=V0X
      VY=V0Y
      VZ=V0Z

cerror 4.7.2018, due to error above, now here
      acc=icharge*echarge1/(gamma*emasskg1)
      vxp=acc*(vy*bzin-vz*byin)
      vyp=acc*(vz*bxin-vx*bzin)
      vzp=acc*(vx*byin-vy*bxin)

999   CONTINUE

      IF (kefield.NE.0) THEN
        dgammao=dgamma
        CALL estep(X,Y,Z,VX,VY,VZ,deltat,gamma,dgamma)
        dgamma=dgammao+dgamma
        VXP=(VX-V0X)/DELTAT
        VYP=(VY-V0Y)/DELTAT
        VZP=(VZ-V0Z)/DELTAT
      ENDIF

      IF (IUSTEP.NE.0) THEN
        dgammao=dgamma
        CALL USTEP(X,Y,Z,VX,VY,VZ,vxp,vyp,vzp,deltat,gamma,dgamma)
        dgamma=dgammao+dgamma
      ENDIF

      RETURN
      End
