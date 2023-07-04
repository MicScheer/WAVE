*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  1.01/00 27/11/97  16.20.29  by  Michael Scheer
*-- Author :    Michael Scheer   27/11/97

      SUBROUTINE BCHARGE(X,Y,Z,BX,BY,BZ,AX,AY,AZ)

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

      IMPLICIT NONE

      INTEGER NCHARP,LUNGFO
      PARAMETER (NCHARP=10)

      INTEGER LUNCHAR,NREAD,ICAL,I

      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,AX,AY,AZ
      DOUBLE PRECISION X0(NCHARP),Y0(NCHARP),Z0(NCHARP),Q(NCHARP),QSUM
      DOUBLE PRECISION XX0,YY0,ZZ0,QQ,RX,RY,RZ,R,R31

      CHARACTER(80) FILE

      DATA LUNCHAR/99/,LUNGFO/16/,ICAL/0/
      DATA FILE/'WI:BCHARGE.DAT'/

      IF (ICAL.EQ.0) THEN

      NREAD=0
      QSUM=0.D0
      OPEN(UNIT=LUNCHAR,FILE=FILE,STATUS='OLD')

11        READ(LUNCHAR,*,END=99)XX0,YY0,ZZ0,QQ
          NREAD=NREAD+1
          X0(NREAD)=XX0
          Y0(NREAD)=YY0
          Z0(NREAD)=ZZ0
          Q(NREAD)=QQ

          QSUM=QSUM+QQ

          GOTO 11

99    CLOSE(LUNCHAR)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SUBROUTINE BCHARGE:'
      WRITE(LUNGFO,*)

      DO I=1,NREAD
         WRITE(LUNGFO,*)'     ',
     &          SNGL(X0(I)),SNGL(Y0(I)),SNGL(Z0(I)),SNGL(Q(I))
      ENDDO
      WRITE(LUNGFO,*)'     QSUM:',QSUM
      WRITE(LUNGFO,*)

      ICAL=1

      ENDIF !ICAL

      BX=0.D0
      BY=0.D0
      BZ=0.D0
      AX=0.D0
      AY=0.D0
      AZ=0.D0
      DO I=1,NREAD

          RX=X-X0(I)
          RY=Y-Y0(I)
          RZ=Z-Z0(I)

          R=DSQRT(RX*RX+RY*RY+RZ*RZ)
          IF (R.EQ.0.D0) THEN
         WRITE(LUNGFO,*)'*** ERROR IN BCHARGE:'
         WRITE(LUNGFO,*)'X,Y,Z COINCIDES WITH POINT CHARGE',I
         WRITE(LUNGFO,*)'X,Y,Z:',X,Y,Z
         WRITE(6,*)'*** ERROR IN BCHARGE:'
         WRITE(6,*)'X,Y,Z COINCIDES WITH POINT CHARGE',I
         WRITE(6,*)'X,Y,Z:',X,Y,Z
         STOP
          ENDIF
          R31=1.D0/(R*R*R)

          BX=BX+Q(I)*RX*R31
          BY=BY+Q(I)*RY*R31
          BZ=BZ+Q(I)*RZ*R31


      ENDDO !NREAD

      RETURN
      END
