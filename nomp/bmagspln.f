*CMZ :  4.00/17 04/10/2022  08.10.22  by  Michael Scheer
*CMZ :  4.00/16 09/09/2022  17.24.46  by  Michael Scheer
*CMZ :  3.02/03 04/11/2014  12.27.16  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.63/03 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.62/02 16/07/2007  09.37.01  by  Michael Scheer
*CMZ :  2.52/16 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.48/04 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.41/13 22/08/2002  13.57.52  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.20/01 08/11/2000  16.17.05  by  Michael Scheer
*CMZ :  2.16/08 31/10/2000  14.25.16  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.25.04  by  Michael Scheer
*CMZ :  2.12/00 03/06/99  10.48.13  by  Michael Scheer
*CMZ :  2.11/00 12/05/99  12.07.51  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.05.53  by  Michael Scheer
*CMZ :  1.00/00 31/07/97  10.47.01  by  Michael Scheer
*CMZ : 00.01/09 31/08/95  15.27.27  by  Michael Scheer
*CMZ : 00.01/08 22/06/95  13.10.15  by  Michael Scheer
*-- Author : Michael Scheer   22/06/95

      SUBROUTINE BMAGSPLN(X,Y,Z,BX,BY,BZ)
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

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEEP,workf90u.
      include 'workf90u.cmn'
*KEND.

C--- USES STORED MAGNETIC FIELD OF REFERENCE ORBIT TO EVALUATE MAGNETIC
C    FIELD BY SPLINE INTERPOLATION

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,debugwave.
      include 'debugwave.cmn'
*KEND.

      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,DUM
      DOUBLE PRECISION XSTOPR
      DOUBLE PRECISION H,A,B,A3A,B3B,H26,H1

      INTEGER ICAL,I,MCODE,MDIM
      INTEGER KLO,KHI,K,KD

      DATA KLO/1/
      DATA ICAL/0/

      IF (ICAL.EQ.0) THEN

        DUM=Y
        DUM=Z

        MDIM=NCO

        IF (IMAGSPLN.EQ.-999) THEN

          IF (IXAMAG_I.EQ.0) THEN
            ALLOCATE(XAMAG(NCO))
            ALLOCATE(BXAMAG(NCO))
            ALLOCATE(BYAMAG(NCO))
            ALLOCATE(BZAMAG(NCO))
            IXAMAG_I=1
          ENDIF

        ELSE IF (IMAGSPLN.LT.0) THEN

          IF (IXAMAG_I.EQ.0) THEN
             ALLOCATE(XAMAG(NCO))
            ALLOCATE(BXAMAG(NCO))
            ALLOCATE(BYAMAG(NCO))
            ALLOCATE(BZAMAG(NCO))
          IXAMAG_I=1
        ENDIF

        OPEN(UNIT=99,FILE='magjob.dat',STATUS='NEW'
     &          ,RECL=256)
        WRITE(99,*)ICODE
          WRITE(99,*)NCO,IBYONLY
          WRITE(99,*)WTRA(1,1,1),WTRA(2,1,1),WTRA(3,1,1)
          WRITE(99,*)WTRA(1,2,1),WTRA(2,2,1),WTRA(3,2,1)
          WRITE(99,*)WTRA(1,1,MDIM),DMYENERGY
          DO I=1,MDIM
            WRITE(99,*)XAMAG(I),BXAMAG(I),BYAMAG(I),BZAMAG(I)
          ENDDO
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '     SR BMAGSPLN:'
          WRITE(LUNGFO,*)
     &      '     MAGNETIC FIELD ARRAY WRITTEN FILE magjob.dat'
          WRITE(LUNGFO,*)
          CLOSE(99)

        ELSE    !IMAGSPLN

          OPEN(UNIT=99,FILE='magjob.dat',STATUS='OLD'
     &      ,RECL=256)

          READ(99,*)MCODE
          MCODE=ABS(MCODE)

          READ(99,*)MDIM,IBYONLY

          IF (IXAMAG_I.EQ.0) THEN
            ALLOCATE(XAMAG(MDIM))
            ALLOCATE(BXAMAG(MDIM))
            ALLOCATE(BYAMAG(MDIM))
            ALLOCATE(BZAMAG(MDIM))
            ALLOCATE(BY2A(MDIM))

            IF (IBYONLY.EQ.0) THEN
              ALLOCATE(BX2A(MDIM))
              ALLOCATE(BZ2A(MDIM))
            ENDIF
            IXAMAG_I=1
          ENDIF

          READ(99,*)XSTART,YSTART,ZSTART
          READ(99,*)VXIN,VYIN,VZIN
          READ(99,*)XSTOPR,DMYENERGY
          IF (XSTOP.EQ.9999.) XSTOP=XSTOPR

          DO I=1,MDIM
            READ(99,*)XAMAG(I),BXAMAG(I),BYAMAG(I),BZAMAG(I)
          ENDDO

          CLOSE(99)

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     BMAGSPLN:'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '     MAGNETIC FIELD READ FROM FILE magjob.dat'
          WRITE(LUNGFO,*)
     &      '     INITIALIZATION OF ELECTRON OVERWRITTEN'
          WRITE(LUNGFO,*)

        ENDIF   !IMAGSPLN

        ALLOCATE(WS1(MDIM))
        ALLOCATE(WS2(MDIM))
        ALLOCATE(WS3(MDIM))
        ALLOCATE(WS4(MDIM))

        IF (IBYONLY.EQ.0) THEN

          CALL util_spline_coef(XAMAG,BXAMAG,MDIM,-9999.0d0,-9999.0d0,BX2A,WS1,WS2,WS3,WS4)
          CALL util_spline_coef(XAMAG,BYAMAG,MDIM,-9999.0d0,-9999.0d0,BY2A,WS1,WS2,WS3,WS4)
          CALL util_spline_coef(XAMAG,BZAMAG,MDIM,-9999.0d0,-9999.0d0,BZ2A,WS1,WS2,WS3,WS4)

        ELSE

          CALL util_spline_coef(XAMAG,BYAMAG,MDIM,-9999.0d0,-9999.0d0,BY2A,WS1,WS2,WS3,WS4)

        ENDIF

        DEALLOCATE(WS1)
        DEALLOCATE(WS2)
        DEALLOCATE(WS3)
        DEALLOCATE(WS4)

        ICAL=1

      ENDIF !ICAL

      IF (IBYONLY.EQ.0) THEN

        IF (X.LT.XAMAG(1)) THEN
          IF (XAMAG(1)-X.LT.2.*(XAMAG(2)-XAMAG(1))) THEN
            BX=BXAMAG(1)+(BXAMAG(2)-BXAMAG(1))/(XAMAG(2)-XAMAG(1))*(X-XAMAG(1))
            BY=BYAMAG(1)+(BYAMAG(2)-BYAMAG(1))/(XAMAG(2)-XAMAG(1))*(X-XAMAG(1))
            BZ=BZAMAG(1)+(BZAMAG(2)-BZAMAG(1))/(XAMAG(2)-XAMAG(1))*(X-XAMAG(1))
            RETURN
          ELSE
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN BMAGSPLN: X OUT OF RANGE ***'
            WRITE(LUNGFO,*)'TRY TO INCREASE NPLOI OR TO DECREASE MYINUM'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN BMAGSPLN: X OUT OF RANGE ***'
            WRITE(6,*)'TRY TO INCREASE NPLOI OR TO DECREASE MYINUM'
            WRITE(6,*)
            STOP
          ENDIF
        ENDIF

        IF (X.GT.XAMAG(MDIM)) THEN
          IF (X-XAMAG(MDIM).LT.2.*(XAMAG(MDIM)-XAMAG(MDIM-1))) THEN
            BX=BXAMAG(MDIM-1)+(BXAMAG(MDIM)-BXAMAG(MDIM-1))/(XAMAG(MDIM)-XAMAG(MDIM-1))*(X-XAMAG(MDIM-1))
            BY=BYAMAG(MDIM-1)+(BYAMAG(MDIM)-BYAMAG(MDIM-1))/(XAMAG(MDIM)-XAMAG(MDIM-1))*(X-XAMAG(MDIM-1))
            BZ=BZAMAG(MDIM-1)+(BZAMAG(MDIM)-BZAMAG(MDIM-1))/(XAMAG(MDIM)-XAMAG(MDIM-1))*(X-XAMAG(MDIM-1))
            RETURN
          ELSE
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN BMAGSPLN: X OUT OF RANGE ***'
            WRITE(LUNGFO,*)'TRY TO INCREASE NPLOI OR TO DECREASE MYINUM'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN BMAGSPLN: X OUT OF RANGE ***'
            WRITE(6,*)'TRY TO INCREASE NPLOI OR TO DECREASE MYINUM'
            WRITE(6,*)
            STOP
          ENDIF
        ENDIF

C--- BMAG_SPLINE_INTER_XYZ{

        IF(     XAMAG(1).LT.XAMAG(MDIM).AND.(X.LT.XAMAG(1).OR.X.GT.XAMAG(MDIM))
     &      .OR.
     &      XAMAG(MDIM).LT.XAMAG(1).AND.(X.LT.XAMAG(MDIM).OR.X.GT.XAMAG(1))) THEN
          STOP '*** ERROR IN BMAGSPLN: X OUT OF RANGE ***'
        ENDIF

        IF (X.GE.XAMAG(KLO)) THEN

C HUNT UP
          KD=1
11        KHI=MIN(KLO+KD,MDIM)
          IF (X.GT.XAMAG(KHI)) THEN
            KD=2*KD
            KLO=KHI
            GOTO 11
          ENDIF

        ELSE    !(X.GE.XAMAG(KLO))

C HUNT DOWN
          KD=1
          KHI=KLO
22        KLO=MAX(KHI-KD,1)
          IF (X.LT.XAMAG(KLO)) THEN
            KD=2*KD
            KHI=KLO
            GOTO 22
          ENDIF

        ENDIF

1       IF (KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XAMAG(K).GT.X)THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
          GOTO 1
        ENDIF

        H=XAMAG(KHI)-XAMAG(KLO)


        IF (H.EQ.0.) THEN
          WRITE(6,*) '*** ERROR IN BMAG_SPLINE_INTER: BAD INPUT ***'
          STOP
        ENDIF

        H1=1.D0/H
        H26=H*H/6.D0
        A=(XAMAG(KHI)-X)*H1
        B=(X-XAMAG(KLO))*H1
        A3A=A*A*A-A
        B3B=B*B*B-B

        BX=A*BXAMAG(KLO)+B*BXAMAG(KHI)+(A3A*BX2A(KLO)+B3B*BX2A(KHI))*H26
        BY=A*BYAMAG(KLO)+B*BYAMAG(KHI)+(A3A*BY2A(KLO)+B3B*BY2A(KHI))*H26
        BZ=A*BZAMAG(KLO)+B*BZAMAG(KHI)+(A3A*BZ2A(KLO)+B3B*BZ2A(KHI))*H26

C--- BMAG_SPLINE_INTER_XYZ}

      ELSE  !IBYONLY

        BX=0.0D0
        BZ=0.0D0

        IF (X.LT.XAMAG(1)) THEN
          IF (XAMAG(1)-X.LT.2.*(XAMAG(2)-XAMAG(1))) THEN
            BY=BYAMAG(1)+(BYAMAG(2)-BYAMAG(1))/(XAMAG(2)-XAMAG(1))*(X-XAMAG(1))
            RETURN
          ELSE
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN BMAGSPLN: X OUT OF RANGE ***'
            WRITE(LUNGFO,*)'TRY TO INCREASE NPLOI OR TO DECREASE MYINUM'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN BMAGSPLN: X OUT OF RANGE ***'
            WRITE(6,*)'TRY TO INCREASE NPLOI OR TO DECREASE MYINUM'
            WRITE(6,*)
            STOP
          ENDIF
        ENDIF

        IF (X.GT.XAMAG(MDIM)) THEN
          IF (X-XAMAG(MDIM).LT.2.*(XAMAG(MDIM)-XAMAG(MDIM-1))) THEN
            BY=BYAMAG(MDIM-1)+(BYAMAG(MDIM)-BYAMAG(MDIM-1))/(XAMAG(MDIM)-XAMAG(MDIM-1))*(X-XAMAG(MDIM-1))
            RETURN
          ELSE
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN BMAGSPLN: X OUT OF RANGE ***'
            WRITE(LUNGFO,*)'TRY TO INCREASE NPLOI OR TO DECREASE MYINUM'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN BMAGSPLN: X OUT OF RANGE ***'
            WRITE(6,*)'TRY TO INCREASE NPLOI OR TO DECREASE MYINUM'
            WRITE(6,*)
            STOP
          ENDIF
        ENDIF

C--- BMAG_SPLINE_INTER_XYZ{

        IF(     XAMAG(1).LT.XAMAG(MDIM).AND.(X.LT.XAMAG(1).OR.X.GT.XAMAG(MDIM))
     &      .OR.
     &      XAMAG(MDIM).LT.XAMAG(1).AND.(X.LT.XAMAG(MDIM).OR.X.GT.XAMAG(1))) THEN
          STOP '*** ERROR IN BMAGSPLN: X OUT OF RANGE ***'
        ENDIF

        IF (X.GE.XAMAG(KLO)) THEN

C HUNT UP
          KD=1
311       KHI=MIN(KLO+KD,MDIM)
          IF (X.GT.XAMAG(KHI)) THEN
            KD=2*KD
            KLO=KHI
            GOTO 311
          ENDIF

        ELSE    !(X.GE.XAMAG(KLO))

C HUNT DOWN
          KD=1
          KHI=KLO
322       KLO=MAX(KHI-KD,1)
          IF (X.LT.XAMAG(KLO)) THEN
            KD=2*KD
            KHI=KLO
            GOTO 322
          ENDIF

        ENDIF

31      IF (KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XAMAG(K).GT.X)THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
          GOTO 31
        ENDIF

        H=XAMAG(KHI)-XAMAG(KLO)

        IF (H.EQ.0.0D0) THEN
          WRITE(6,*) '*** ERROR IN BMAG_SPLINE_INTER: BAD INPUT ***'
          STOP
        ENDIF

        H1=1.0D0/H
        H26=H*H/6.D0
        A=(XAMAG(KHI)-X)*H1
        B=(X-XAMAG(KLO))*H1
        A3A=A*A*A-A
        B3B=B*B*B-B

        BY=A*BYAMAG(KLO)+B*BYAMAG(KHI)+(A3A*BY2A(KLO)+B3B*BY2A(KHI))*H26

C--- BMAG_SPLINE_INTER_XYZ}

      ENDIF !IBYONLY

      RETURN
      END
