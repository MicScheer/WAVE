*CMZ :  2.16/04 18/09/2013  12.33.23  by  Michael Scheer
*CMZ : 00.02/00 07/11/96  17.02.37  by  Michael Scheer
*CMZ : 00.01/09 25/10/95  17.37.17  by  Michael Scheer
*-- Author :    Michael Scheer   29/09/95

      SUBROUTINE BPOLY3D(XIN,YIN,ZIN,BX,BY,BZ,AX,AY,AZ)
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

C--- ALL INDICES ACCORDING TO FORTRAN, BUT LORD3D AND MORD3D REFER TO MATH. INDICES

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,bpoly3d.
      include 'bpoly3d.cmn'
*KEND.

      INTEGER IX,IY,IZ,NPOWP,ICAL
     &         ,IWARNXMN,IWARNXMX,IWARNYMN,IWARNYMX,IWARNZMN,IWARNZMX
     &         ,IBX,IBY,IBZ
     &         ,IND

      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,XPOW1,YPOW1,ZPOW1,AX,AY,AZ
     &                  ,XIN,YIN,ZIN,DLEN,DLEN2

      PARAMETER (NPOWP=NDIMC+1)
      DOUBLE PRECISION XPOW(NPOWP),YPOW(NPOWP),ZPOW(NPOWP)

      DATA ICAL/0/
      DATA IWARNXMN,IWARNXMX,IWARNYMN,IWARNYMX,IWARNZMN,IWARNZMX/6*0/

      IF (ICAL.EQ.0) THEN
          CALL BPOLY3DINI
          DLEN2=2.D0*X3DMAX
          DLEN=2.D0*DLEN2
          ICAL=1
      ENDIF

      IF (MORD3D.GT.NDIMC-1) THEN
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BPOLY3D: DIMENSION NDIMC EXCEEDED ***'
          WRITE(6,*)
          STOP
      ENDIF

      IF (MORD3D.GT.NPOWP) THEN
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BPOLY3D: DIMENSION NPOWP EXCEEDED ***'
          WRITE(6,*)
          STOP
      ENDIF

      X3DMIN=X3DMIN*B3DSCALE
      X3DMAX=X3DMAX*B3DSCALE
      Y3DMIN=Y3DMIN*B3DSCALE
      Y3DMAX=Y3DMAX*B3DSCALE
      Z3DMIN=Z3DMIN*B3DSCALE
      Z3DMAX=Z3DMAX*B3DSCALE

      X=XIN*B3DSCALE
      Y=YIN*B3DSCALE
      Z=ZIN*B3DSCALE

      IBX=1
      IBY=1
      IBZ=1
      IF (I3DQUART.NE.0) THEN

          X=DMOD(X,DLEN)

          IF (X.GT.DLEN2) THEN
         X=X-DLEN
          ELSEIF (X.LT.-DLEN2) THEN
         X=X+DLEN
          ENDIF

          IF(X.LT.0.0D0) THEN
         X=-X
         IBX=-IBX
          ENDIF

          IF(Y.LT.0.0D0) THEN
         Y=-Y
         IBZ=-IBZ
          ENDIF

          IF(Z.LT.0.0) THEN
         Z=-Z
         IBZ=-IBZ
          ENDIF

          IF(X.GT.X3DMAX) THEN
         X=X3DMAX-(X-X3DMAX)
         IBY=-IBY
         IBZ=-IBZ
          ENDIF

          IF(X.LT.0.0D0) THEN
         X=-X
         IBX=-IBX
          ENDIF

      ENDIF !I3DQUART

      IF (IWARNXMN.EQ.0.AND.X.LT.X3DMIN) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR BPOLY3D: X OUT OF FIT RANGE'
          WRITE(LUNGFO,*)'X .LT. X3DMIN'
          WRITE(LUNGFO,*)'X, X3DMIN:',X,X3DMIN
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR BPOLY3D: X OUT OF FIT RANGE'
          WRITE(6,*)'X .LT. X3DMIN'
          WRITE(6,*)'X, X3DMIN:',X,X3DMIN
          WRITE(6,*)
          IWARNXMN=1
      ENDIF

      IF (IWARNXMX.EQ.0.AND.X.GT.X3DMAX) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR BPOLY3D: X OUT OF FIT RANGE'
          WRITE(LUNGFO,*)'X .GT. X3DMAX'
          WRITE(LUNGFO,*)'X, X3DMAX:',X,X3DMAX
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR BPOLY3D: X OUT OF FIT RANGE'
          WRITE(6,*)'X .GT. X3DMAX'
          WRITE(6,*)'X, X3DMAX:',X,X3DMAX
          WRITE(6,*)
          IWARNXMX=1
      ENDIF

      IF (IWARNYMN.EQ.0.AND.Y.LT.Y3DMIN) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR BPOLY3D: Y OUT OF FIT RANGE'
          WRITE(LUNGFO,*)'Y .LT. Y3DMIN'
          WRITE(LUNGFO,*)'Y, Y3DMIN:',Y,Y3DMIN
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR BPOLY3D: Y OUT OF FIT RANGE'
          WRITE(6,*)'Y .LT. Y3DMIN'
          WRITE(6,*)'Y, Y3DMIN:',Y,Y3DMIN
          WRITE(6,*)
          IWARNYMN=1
      ENDIF

      IF (IWARNYMX.EQ.0.AND.Y.GT.Y3DMAX) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR BPOLY3D: Y OUT OF FIT RANGE'
          WRITE(LUNGFO,*)'Y .GT. Y3DMAX'
          WRITE(LUNGFO,*)'Y, Y3DMAX:',Y,Y3DMAX
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR BPOLY3D: Y OUT OF FIT RANGE'
          WRITE(6,*)'Y .GT. Y3DMAX'
          WRITE(6,*)'Y, Y3DMAX:',Y,Y3DMAX
          WRITE(6,*)
          IWARNYMX=1
      ENDIF

      IF (IWARNZMN.EQ.0.AND.Z.LT.Z3DMIN) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR BPOLY3D: Z OUT OF FIT RANGE'
          WRITE(LUNGFO,*)'Z .LT. Z3DMIN'
          WRITE(LUNGFO,*)'Z, Z3DMIN:',Z,Z3DMIN
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR BPOLY3D: Z OUT OF FIT RANGE'
          WRITE(6,*)'Z .LT. Z3DMIN'
          WRITE(6,*)'Z, Z3DMIN:',Z,Z3DMIN
          WRITE(6,*)
          IWARNZMN=1
      ENDIF

      IF (IWARNZMX.EQ.0.AND.Z.GT.Z3DMAX) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR BPOLY3D: Z OUT OF FIT RANGE'
          WRITE(LUNGFO,*)'Z .GT. Z3DMAX'
          WRITE(LUNGFO,*)'Z, Z3DMAX:',Z,Z3DMAX
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR BPOLY3D: Z OUT OF FIT RANGE'
          WRITE(6,*)'Z .GT. Z3DMAX'
          WRITE(6,*)'Z, Z3DMAX:',Z,Z3DMAX
          WRITE(6,*)
          IWARNZMX=1
      ENDIF

      XPOW(1)=1.D0
      DO IX=2,MORD3D
          XPOW(IX)=XPOW(IX-1)*X
      ENDDO

      YPOW(1)=1.D0
      DO IY=2,MORD3D
          YPOW(IY)=YPOW(IY-1)*Y
      ENDDO

      ZPOW(1)=1.D0
      DO IZ=2,MORD3D
          ZPOW(IZ)=ZPOW(IZ-1)*Z
      ENDDO

      AX=0.D0
      AY=0.D0
      AZ=0.D0

      BX=0.D0
      BY=0.D0
      BZ=0.D0

      DO IND=1,NCIND

             IX=ICIND(1,IND)
             IY=ICIND(2,IND)
             IZ=ICIND(3,IND)

             IF (IX.GT.1) THEN
            XPOW1=XPOW(IX-1)
             ELSE
            XPOW1=0.D0
             ENDIF
             BX=BX-(IX-1)*CIND(IND)*XPOW1*YPOW(IY)*ZPOW(IZ)

             IF (IY.GT.1) THEN
            YPOW1=YPOW(IY-1)
             ELSE
            YPOW1=0.D0
             ENDIF

             BY=BY-(IY-1)*CIND(IND)*XPOW(IX)*YPOW1*ZPOW(IZ)

             IF (IZ.GT.1) THEN
            ZPOW1=ZPOW(IZ-1)
             ELSE
            ZPOW1=0.D0
             ENDIF

             BZ=BZ-(IZ-1)*CIND(IND)*XPOW(IX)*YPOW(IY)*ZPOW1

      ENDDO

      IF (IBX.LT.0) BX=-BX
      IF (IBY.LT.0) BY=-BY
      IF (IBZ.LT.0) BZ=-BZ

      X3DMIN=X3DMIN/B3DSCALE
      X3DMAX=X3DMAX/B3DSCALE
      Y3DMIN=Y3DMIN/B3DSCALE
      Y3DMAX=Y3DMAX/B3DSCALE
      Z3DMIN=Z3DMIN/B3DSCALE
      Z3DMAX=Z3DMAX/B3DSCALE

      RETURN
      END

