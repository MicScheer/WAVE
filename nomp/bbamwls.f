*CMZ :  4.00/04 17/05/2019  14.17.20  by  Michael Scheer
*CMZ :  2.63/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.56/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.56/01 21/10/2005  12.22.07  by  Michael Scheer
*CMZ :  2.54/05 19/05/2005  17.01.55  by  Michael Scheer
*-- Author :    Michael Scheer   22/11/96
      SUBROUTINE BBAMWLS(X,Y,Z,BX,BY,BZ)
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

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,bamwls.
      include 'bamwls.cmn'
*KEEP,fourier.
      include 'fourier.cmn'
*KEND.

      DOUBLE PRECISION
     &  X,Y,Z,BX,BY,BZ,BXX,BYY,BZZ,
     &  XK0FOUR1,ZK0FOUR1,ZL0FOUR1,A01,A1(MAXFOUR),
     &  XKFOUR1(MAXFOUR),YKFOUR1(MAXFOUR),ZKFOUR1(MAXFOUR),
     &  XK0FOUR2,ZK0FOUR2,ZL0FOUR2,A02,A2(MAXFOUR),
     &  XKFOUR2(MAXFOUR),YKFOUR2(MAXFOUR),ZKFOUR2(MAXFOUR),
     &  AXX,AYY,AZZ,B01,B02,
     &  P1(4),B1,B2

      INTEGER ICAL,I,IK,IRFILB0O,LUNDUM,NFOUR1,NFOUR2

      DATA ICAL/0/,LUNDUM/99/
      DATA P1/-0.00150191,0.00508334,7.2303E-07,-2.2247E-08/

      IF (ICAL.EQ.0) THEN

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      BBAMWLS:'
        WRITE(LUNGFO,*)'      MODE:',MBAMWLS
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      CORRMS,CORRMM:',SNGL(CORRMS),SNGL(CORRMM)
        WRITE(LUNGFO,*)'      (scaling factors for neg. and. pos. field of WLS)'
        WRITE(LUNGFO,*)

C ERSTER STEERER

        OPEN(UNIT=LUNDUM,FILE='bamwls_neuer_steerer.fou',status='old')

        READ(LUNDUM,*)
        READ(LUNDUM,*)ZL0FOUR1
        READ(LUNDUM,*)NFOUR1

        IF (NFOUR1.GT.MAXFOUR) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '*** ERROR IN BFOUR ***'
          WRITE(LUNGFO,*)'NFOUR1.GT.MAXFOUR'
          WRITE(6,*)
          WRITE(6,*)
     &      '*** ERROR IN BFOUR ***'
          WRITE(6,*)'NFOUR1.GT.MAXFOUR'
          STOP
        ENDIF

        READ(LUNDUM,*)IK,A01

        DO I=1,NFOUR1-1
          READ(LUNDUM,*)IK,A1(IK-1)
        END DO

        CLOSE(LUNDUM)

        IF (XLCORRL.NE.0.0D0) THEN
          XK0FOUR1=2.D0*PI1/XLCORRL
        ELSE
          XK0FOUR1=0.0D0
        ENDIF

        ZK0FOUR1=2.D0*PI1/ZL0FOUR1

        B01=A01/2.0D0

        DO I=1,NFOUR1-1
          ZKFOUR1(I)=ZK0FOUR1*I
          XKFOUR1(I)=XK0FOUR1 !VORERST, SIEHE AUCH OBEN
          YKFOUR1(I)=DSQRT(ZKFOUR1(I)**2+XKFOUR1(I)**2)
          B01=B01+A1(I)
        END DO

        B1=P1(1)+CURRL*(P1(2)+CURRL*(P1(3)+P1(4)*CURRL))
        A01=A01/B01*B1*CORRL
        A1=A1/B01*B1*CORRL

C ZWEITER STEERER

        OPEN(UNIT=LUNDUM,FILE='bamwls_alter_steerer.fou',status='old')

        READ(LUNDUM,*)
        READ(LUNDUM,*)ZL0FOUR2
        READ(LUNDUM,*)NFOUR2

        IF (NFOUR2.GT.MAXFOUR) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '*** ERROR IN BFOUR ***'
          WRITE(LUNGFO,*)'NFOUR2.GT.MAXFOUR'
          WRITE(6,*)
          WRITE(6,*)
     &      '*** ERROR IN BFOUR ***'
          WRITE(6,*)'NFOUR2.GT.MAXFOUR'
          STOP
        ENDIF

        READ(LUNDUM,*)IK,A02

        DO I=1,NFOUR2-1
          READ(LUNDUM,*)IK,A2(IK-1)
        END DO

        CLOSE(LUNDUM)

        IF (XLCORRR.NE.0.0D0) THEN
          XK0FOUR2=2.D0*PI1/XLCORRR
        ELSE
          XK0FOUR2=0.0D0
        ENDIF

        ZK0FOUR2=2.D0*PI1/ZL0FOUR2

        DO I=1,NFOUR2-1
          ZKFOUR2(I)=ZK0FOUR2*I
          XKFOUR2(I)=XK0FOUR2 !VORERST, SIEHE AUCH OBEN
          YKFOUR2(I)=DSQRT(ZKFOUR2(I)**2+XKFOUR2(I)**2)
          B02=B02+A2(I)
        END DO

        B2=(CURRR*(-1.0957*CURRR+215.31)-6.05)/10000.
        B02=(14.436*(-1.0957*14.436+215.31)-6.05)/10000.
        A02=A02/B02*B2*CORRR
        A2=A2/B02*B2*CORRR

        WRITE(LUNGFO,*)'      CURRL,CURRR:    ',SNGL(CURRL),SNGL(CURRR)
        WRITE(LUNGFO,*)'      CORRL,CORRR:    ',SNGL(CORRL),SNGL(CORRR)
        WRITE(LUNGFO,*)'      max. fields:    ',SNGL(B1),SNGL(B2)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      XCORRL,XCORRR:  ',SNGL(XCORRL),SNGL(XCORRR)
        WRITE(LUNGFO,*)'      XLCORRL,XLCORRR:',SNGL(XLCORRL),SNGL(XLCORRR)
        WRITE(LUNGFO,*)

        ICAL=1

      ENDIF

      IRFILB0O=IRFILB0

      BX=0.0D0
      BY=0.0D0
      BZ=0.0D0

C STEERER

      CALL BFOURMULT(X-XCORRL,Y,Z,BXX,BYY,BZZ,
     &  AXX,AYY,AZZ,
     &  NFOUR1,ZL0FOUR1,XKFOUR1,YKFOUR1,ZKFOUR1,A01,A1)

      BX=BX+BXX
      BY=BY+BYY
      BZ=BZ+BZZ

      CALL BFOURMULT(X-XCORRR,Y,Z,BXX,BYY,BZZ,
     &  AXX,AYY,AZZ,
     &  NFOUR2,ZL0FOUR1,XKFOUR2,YKFOUR2,ZKFOUR2,A02,A2)

      BX=BX+BXX
      BY=BY+BYY
      BZ=BZ+BZZ

      IF (MBAMWLS.EQ.0) THEN

C RADIA-FELDMAPPE DES WLS

        IRFILB0=-2
        CALL BMESS(ABS(X),ABS(Y),ABS(Z),BXX,BYY,BZZ)

        IF (X.LT.0.0D0) THEN
          BXX=-BXX
        ENDIF
        IF (Y.LT.0.0D0) THEN
          BXX=-BXX
          BZZ=-BZZ
        ENDIF
        IF (Z.LT.0.0D0) THEN
          BZZ=-BZZ
        ENDIF

      ELSE IF (MBAMWLS.EQ.1) THEN

        CALL BTAB(X,Y,Z,BXX,BYY,BZZ,AXX,AYY,AZZ)

      ENDIF !(MBAMWLS.EQ.0) THEN

      IF (BYY.LT.0.0D0) THEN
        BXX=BXX*CORRMS
        BYY=BYY*CORRMS
        BZZ=BZZ*CORRMS
      ELSE
        BXX=BXX*CORRMM
        BYY=BYY*CORRMM
        BZZ=BZZ*CORRMM
      ENDIF

      BX=BX+BXX
      BY=BY+BYY
      BZ=BZ+BZZ

      IRFILB0=IRFILB0O

      RETURN
      END
