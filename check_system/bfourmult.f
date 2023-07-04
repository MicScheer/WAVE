*CMZ :  2.55/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.54/05 19/05/2005  14.58.44  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BFOURMULT(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,
     &  AXOUT,AYOUT,AZOUT,
     &  NFOUR,ZL0FOUR,XKFOUR,YKFOUR,ZKFOUR,A0,A)
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

C SEE BFOUR

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      COMPLEX*16 CDEXPOMX,CEXPOMZ,CDEXPOMZ

      DOUBLE PRECISION A0,A(MAXFOUR),
     &  XKFOUR(MAXFOUR),YKFOUR(MAXFOUR),ZKFOUR(MAXFOUR),
     &  ZL0FOUR

      DOUBLE PRECISION XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT,
     &  DSNXKX,DCSXKX,DSHYKY,DCHYKY,DSNZKZ,DCSZKZ
     &  ,BXH,BYH,BZH,AXH,AYH,AZH,X

      DOUBLE PRECISION EXPOMY,DEXPOMY,EXPOMY1

      DOUBLE PRECISION DNULL

      INTEGER K,NFOUR

      DATA DNULL/0.D0/

      IF (NFOUR.NE.-9999.AND.(NFOUR.LT.1.OR.NFOUR.GT.MAXFOUR)) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN BFOURMULT ***'
        WRITE(LUNGFO,*)'NFOUR.LT.1.OR.NFOUR.GT.MAXFOUR'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN BFOURMULT ***'
        WRITE(6,*)'NFOUR.LT.1.OR.NFOUR.GT.MAXFOUR'
        WRITE(6,*)
        STOP
      ENDIF

      IF (DABS(XIN).GT.ZL0FOUR/2.0D0) THEN
        BXOUT=0.0D0
        BYOUT=0.0D0
        BZOUT=0.0D0
        AXOUT=0.0D0
        AYOUT=0.0D0
        AZOUT=0.0D0
        RETURN
      ELSE
        X=XIN
      ENDIF

      BXH=0.
      BYH=A0/2.D0
      BZH=0.

C IF CHANGED, CONSIDER FOLLOWING LOOP AND SR MYBFELD {



      AXH= A0/2.D0*  XIN !XIN IS HERE Z
      AYH=0.
      AZH=0.

C IF CHANGED, CONSIDER FOLLOWING LOOP AND SR MYBFELD }

      CDEXPOMX=CDEXP(DCMPLX(DNULL,XKFOUR(1)*(-ZIN)))
      DCSXKX=DREAL(CDEXPOMX)
      DSNXKX=DIMAG(CDEXPOMX)

      DEXPOMY=DEXP(YKFOUR(1)*YIN)
      EXPOMY=1.D0

      CDEXPOMZ=CDEXP(DCMPLX(DNULL,ZKFOUR(1)*    X ))
      CEXPOMZ=DCMPLX(1.D0,DNULL)

      DO K=1,NFOUR-1

        IF (XKFOUR(1).NE.0.0D0) THEN
          EXPOMY=DEXP(YKFOUR(K)*YIN)
        ELSE
          EXPOMY=EXPOMY*DEXPOMY
        ENDIF
        EXPOMY1=1.D0/EXPOMY
        DCHYKY=(EXPOMY+EXPOMY1)*0.5D0
        DSHYKY=(EXPOMY-EXPOMY1)*0.5D0

        CEXPOMZ=CEXPOMZ*CDEXPOMZ
        DCSZKZ=DREAL(CEXPOMZ)
        DSNZKZ=DIMAG(CEXPOMZ)

        BXH=BXH-A(K)*XKFOUR(K)/YKFOUR(K)*DSNXKX*DSHYKY*DCSZKZ
        BYH=BYH+A(K)*                    DCSXKX*DCHYKY*DCSZKZ
        BZH=BZH-A(K)*ZKFOUR(K)/YKFOUR(K)*DCSXKX*DSHYKY*DSNZKZ



        AXH=AXH+A(K)/ZKFOUR(K)*DCSXKX*DCHYKY*DSNZKZ
        AZH=AZH+0.0

        AYH=AYH+A(K)/ZKFOUR(K)*XKFOUR(K)/YKFOUR(K)*DSNXKX*DSHYKY*DSNZKZ

      ENDDO

      BZOUT=-BXH
      BYOUT= BYH
      BXOUT= BZH

      AZOUT=-AXH
      AYOUT= AYH
      AXOUT= AZH

      RETURN
      END
