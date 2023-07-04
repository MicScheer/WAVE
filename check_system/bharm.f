*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  1.03/06 11/06/98  13.46.47  by  Michael Scheer
*CMZ : 00.01/08 04/04/95  17.31.58  by  Michael Scheer
*CMZ : 00.01/03 28/11/94  13.20.41  by  Michael Scheer
*CMZ :  0.00/03 14/11/94  12.38.33  by  Michael Scheer
*CMZ :  0.00/02 02/11/94  11.28.11  by  Michael Scheer
*CMZ :  0.00/01 31/10/94  09.48.00  by  Michael Scheer
*CMZ :  0.00/00 28/10/94  16.14.51  by  Michael Scheer
*-- Author :    Michael Scheer   28/10/94
C***************************************************************
        SUBROUTINE BHARM(X,Y,Z,BX,BY,BZ
     &  ,NFIRSTX,NORDX,NSTEPX
     &  ,NFIRSTY,NORDY,NSTEPY
     &  ,Q,QA0,QA,NQDIM,XKX,YKY,ZKZ,IFHALBA,GAP2PI,WIDTH)

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

      INTEGER NDUMP
      PARAMETER (NDUMP=1000)

      INTEGER NFIRSTX,NORDX,NSTEPX,NFIRSTY,NORDY,NSTEPY,NQDIM
      INTEGER IORDX,IORDY
      INTEGER IFHALBA

        DOUBLE PRECISION DUM,DXC,DXS

      DOUBLE PRECISION GAP2PI,WIDTH

      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,Q(NQDIM,NQDIM),XKX,YKY,ZKZ
      DOUBLE PRECISION QA0(NQDIM),QA(NQDIM,NQDIM)

      DOUBLE PRECISION XSIN(NDUMP),XCOS(NDUMP)
      DOUBLE PRECISION YSINH(NDUMP),YCOSH(NDUMP)
      DOUBLE PRECISION ZSIN(NDUMP),ZCOS(NDUMP)
      DOUBLE PRECISION ZIJC,ZIJS,YIJC,YIJS,XKIJ,ZKIJ,RKIJ


      IF (NQDIM.GT.NDUMP) THEN
      WRITE(6,*)'*** ERROR IN BHARM: DIMENSION EXCEEDED ***'
      STOP
      ENDIF

      IF (IFHALBA.NE.0) THEN

      CALL B_FIELD(X,Y,Z,Bx,By,Bz,QA0,QA,GAP2PI,ZKZ,width,
     &             NFIRSTX,NORDX,NSTEPX,NFIRSTY,NORDY,NSTEPY,NQDIM)

      ELSE  !IFHALBA

      BX=0.0
      BY=0.0
      BZ=0.0

      IF (YKY.GT.0.0) THEN

      WRITE(6,*)'ADDITIONSTHEOREME EINBAUEN'
        DO IORDY=NFIRSTY,NORDY,NSTEPY
          YCOSH(IORDY)=DCOSH(Y*(IORDY-1)*YKY)
          YSINH(IORDY)=DSINH(Y*(IORDY-1)*YKY)
        ENDDO

        DO IORDX=NFIRSTX,NORDX,NSTEPX
          XCOS(IORDX)=DCOS(X*(IORDX-1)*XKX)
          XSIN(IORDY)=DSIN(X*(IORDX-1)*XKX)
        ENDDO

        DO IORDX=NFIRSTX,NORDX,NSTEPX
        DO IORDY=NFIRSTY,NORDY,NSTEPY

          IF ((((IORDY-1)*YKY)**2-((IORDX-1)*XKX)**2).GE.0.0) THEN
            RKIJ=(IORDY-1)*YKY
            XKIJ=(IORDX-1)*XKX
            ZKIJ=DSQRT(((IORDY-1)*YKY)**2-((IORDX-1)*XKX)**2)
            ZIJC=DCOS (Z*ZKIJ)
            ZIJS=DSIN (Z*ZKIJ)
          ELSE
            RKIJ=1.E-30
            XKIJ=0.0
            ZKIJ=0.0
            YCOSH(IORDY)=0.0
            YSINH(IORDY)=0.0
            ZIJC=0.0
            ZIJS=0.0
          ENDIF

          BX=BX-XKIJ/RKIJ*Q(IORDX,IORDY)*XSIN(IORDX)*YSINH(IORDY)*ZIJC
          BY=BY+          Q(IORDX,IORDY)*XCOS(IORDX)*YCOSH(IORDY)*ZIJC
          BZ=BZ-ZKIJ/RKIJ*Q(IORDX,IORDY)*XCOS(IORDX)*YSINH(IORDY)*ZIJS

        ENDDO
        ENDDO

      ENDIF   !YKY>0

      IF (YKY.LT.0.0) THEN


        DUM=X*(NFIRSTX-1)*XKX
        XCOS(NFIRSTX)=DCOS(DUM)
        XSIN(NFIRSTX)=DSIN(DUM)
        DUM=X*NSTEPX*XKX
        DXC=DCOS(DUM)
        DXS=DSIN(DUM)

        DO IORDX=NFIRSTX+NSTEPX,NORDX,NSTEPX
          XCOS(IORDX)=XCOS(IORDX-NSTEPX)*DXC-XSIN(IORDX-NSTEPX)*DXS
          XSIN(IORDX)=XCOS(IORDX-NSTEPX)*DXS+XSIN(IORDX-NSTEPX)*DXC
        ENDDO

C        DO IORDX=NFIRSTX,NORDX,NSTEPX
C          XCOS(IORDX)=DCOS (X*(IORDX-1)*XKX)
C          XSIN(IORDX)=DSIN (X*(IORDX-1)*XKX)
C        ENDDO

        DUM=Z*(NFIRSTY-1)*ZKZ
        ZCOS(NFIRSTY)=DCOS(DUM)
        ZSIN(NFIRSTY)=DSIN(DUM)
        DUM=Z*NSTEPY*ZKZ
        DXC=DCOS(DUM)
        DXS=DSIN(DUM)

        DO IORDY=NFIRSTY+NSTEPY,NORDY,NSTEPY
          ZCOS(IORDY)=ZCOS(IORDY-NSTEPY)*DXC-ZSIN(IORDY-NSTEPY)*DXS
          ZSIN(IORDY)=ZCOS(IORDY-NSTEPY)*DXS+ZSIN(IORDY-NSTEPY)*DXC
        ENDDO

C        DO IORDY=NFIRSTY,NORDY,NSTEPY
C          ZCOS(IORDY)=DCOS (Z*(IORDY-1)*ZKZ)
C          ZSIN(IORDY)=DSIN (Z*(IORDY-1)*ZKZ)
C        ENDDO

        DO IORDX=NFIRSTX,NORDX,NSTEPX
        DO IORDY=NFIRSTY,NORDY,NSTEPY
          RKIJ=DSQRT(((IORDY-1)*ZKZ)**2+((IORDX-1)*XKX)**2)
          XKIJ=(IORDX-1)*XKX
          ZKIJ=(IORDY-1)*ZKZ
          YIJS=DSINH(Y*RKIJ)
C          YIJC=DCOSH(Y*RKIJ)
          YIJC=DSQRT(1.D0+YIJS**2)
          BX=BX-XKIJ/RKIJ*Q(IORDX,IORDY)*XSIN(IORDX)*YIJS*ZCOS(IORDY)
          BY=BY+          Q(IORDX,IORDY)*XCOS(IORDX)*YIJC*ZCOS(IORDY)
          BZ=BZ-ZKIJ/RKIJ*Q(IORDX,IORDY)*XCOS(IORDX)*YIJS*ZSIN(IORDY)
        ENDDO
        ENDDO

      ELSEIF (YKY.EQ.0.0) THEN
      WRITE(6,*)'*** ERROR IN BHARM: YKY=0.0 ?? ***'
        STOP
      ENDIF   !YKY.LT.0

      RETURN
      ENDIF !IFHALBA

      END
