*CMZ :  2.56/00 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.52/16 06/01/2005  15.27.37  by  Michael Scheer
*CMZ :  2.52/14 22/12/2004  16.09.45  by  Michael Scheer
*CMZ :  2.52/13 16/12/2004  21.14.00  by  Michael Scheer
*-- Author :    Michael Scheer   15/12/2004
      SUBROUTINE SECTMAGS
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

C--- Approximate vertical magnetic field as sequence of sector magnets

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,track0.
      include 'track0.cmn'
*KEND.

      DOUBLE PRECISION RHO,SPHI1,PHI1,SPHI2,PHI2,DPHI,
     &  X1,Z1,X2,Z2,ZERO,BY,
     &  HMATT(2,2),VMATT(2,2),VMAT(2,2),HMAT(2,2),HFOC,VFOC,WORK22(2,2)

      INTEGER ISEC,JSEC,NSEC,IX,NSTEP,IZERO

      CHARACTER(4) CNAME,CNAME10(10)
      CHARACTER(2) CTYPE

      DATA ZERO/0.0D0/

      IF (ICHARGE.LE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING in SECTMAGS: ICHARGE NEGATIVE'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING in SECTMAGS: ICHARGE NEGATIVE'
        WRITE(6,*)
      ENDIF

      IF (IWSECTMAGS.GT.999) THEN
        WRITE(LUNGFO,*)'*** Error in SECTMAGS: Too many magnets'
        WRITE(LUNGFO,*)'*** decrease IWSECTMAGS in wave.in'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** Programm WAVE aborted ***'
        WRITE(6,*)'*** Error in SECTMAGS: Too many magnets'
        WRITE(6,*)'*** decrease IWSECTMAGS in wave.in'
        WRITE(6,*)
        WRITE(6,*)'*** Programm WAVE aborted ***'
        STOP
      ENDIF

      NSEC=IWSECTMAGS
      NSTEP=(NCO-1)/IWSECTMAGS

      IX=1

      X1=WSXYZ(1,IX)
      Z1=WSXYZ(3,IX)

      PHI1=ATAN2(WVXYZ(3,IX),WVXYZ(1,IX))
      SPHI1=SIN(PHI1)

      ISEC=0

      HMATT(1,1)=1.0D0
      HMATT(1,2)=0.0D0
      HMATT(2,1)=0.0D0
      HMATT(2,2)=1.0D0

      VMATT=HMATT

      OPEN(UNIT=99,FILE='wave.smag',status='new')
      OPEN(UNIT=98,FILE='wave.xmag',status='new')

      WRITE(99,*)ICODE
      WRITE(99,*)CODE
      WRITE(98,*)ICODE
      WRITE(98,*)CODE

      IZERO=0

1     CONTINUE

      BY=WBXYZ(2,IX)
      IF (ABS(BY).LT.1.0D-6.AND.IX.LT.NCO) THEN
        IX=IX+1
        IZERO=1
        GOTO 1
      ELSE
        IX=IX+NSTEP-IZERO
      ENDIF

        ISEC=ISEC+1

        IF (IX.GE.NCO.OR.ISEC.EQ.IWSECTMAGS) THEN
          IX=NCO
          X2=XF0
          Z2=ZF0
          PHI2=ATAN2(VZF0,VXF0)
        ELSE
          X2=WSXYZ(1,IX)
          Z2=WSXYZ(3,IX)
          PHI2=ATAN2(WVXYZ(3,IX),WVXYZ(1,IX))
        ENDIF

        SPHI2=SIN(PHI2)
        DPHI=PHI2-PHI1

        IF (IZERO.EQ.0.AND.SPHI2.NE.SPHI1) THEN

          RHO=((X2-X1)/(SPHI2-SPHI1))

          CNAME='ASEC'
          IF (ISEC.LE.9) THEN
            WRITE(CNAME(4:4),'(I1)') ISEC
          ELSE IF (ISEC.LE.99) THEN
            WRITE(CNAME(3:4),'(I2)') ISEC
          ELSE IF (ISEC.LE.999) THEN
            WRITE(CNAME(2:4),'(I3)') ISEC
          ELSE
            WRITE(CNAME,'(I4)') ISEC
          ENDIF
          CTYPE='CO'
          WRITE(99,'(1X,A4,1X,A2,2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6)')
     &      CNAME,CTYPE,-PHI1,RHO,ZERO,ZERO
          WRITE(98,'(4E16.6)')X1,X2,RHO,-PHI1

          CALL CORNER_MAT(-PHI1,RHO,HMAT,VMAT,0.0D0)
          CALL UTIL_MATRIX_MULTIPLICATION(2,2,2,HMATT,HMAT,HMATT,WORK22)
          CALL UTIL_MATRIX_MULTIPLICATION(2,2,2,VMATT,VMAT,VMATT,WORK22)

          CNAME='BSEC'
          IF (ISEC.LE.9) THEN
            WRITE(CNAME(4:4),'(I1)') ISEC
          ELSE IF (ISEC.LE.99) THEN
            WRITE(CNAME(3:4),'(I2)') ISEC
          ELSE IF (ISEC.LE.999) THEN
            WRITE(CNAME(2:4),'(I3)') ISEC
          ELSE
            WRITE(CNAME,'(I4)') ISEC
          ENDIF
          CTYPE='DI'
          WRITE(99,'(1X,A4,1X,A2,2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6)')
     &      CNAME,CTYPE,DPHI,RHO,ZERO,ZERO
          WRITE(98,'(4E16.6)')X1,X2,RHO,DPHI

          CALL SECTMAG_MAT(DPHI,RHO,HMAT,VMAT)
          CALL UTIL_MATRIX_MULTIPLICATION(2,2,2,HMATT,HMAT,HMATT,WORK22)
          CALL UTIL_MATRIX_MULTIPLICATION(2,2,2,VMATT,VMAT,VMATT,WORK22)

          CNAME='CSEC'
          IF (ISEC.LE.9) THEN
            WRITE(CNAME(4:4),'(I1)') ISEC
          ELSE IF (ISEC.LE.99) THEN
            WRITE(CNAME(3:4),'(I2)') ISEC
          ELSE IF (ISEC.LE.999) THEN
            WRITE(CNAME(2:4),'(I3)') ISEC
          ELSE
            WRITE(CNAME,'(I4)') ISEC
          ENDIF
          CTYPE='CO'
          WRITE(99,'(1X,A4,1X,A2,2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6)')
     &      CNAME,CTYPE,PHI2,RHO,ZERO,ZERO
          WRITE(98,'(4E16.6)')X1,X2,RHO,PHI2

          CALL CORNER_MAT(PHI2,RHO,HMAT,VMAT,0.0D0)
          CALL UTIL_MATRIX_MULTIPLICATION(2,2,2,HMATT,HMAT,HMATT,WORK22)
          CALL UTIL_MATRIX_MULTIPLICATION(2,2,2,VMATT,VMAT,VMATT,WORK22)

        ELSE !SPHI2.NE.SPHI1

          CNAME='DSEC'
          IF (ISEC.LE.9) THEN
            WRITE(CNAME(4:4),'(I1)') ISEC
          ELSE IF (ISEC.LE.99) THEN
            WRITE(CNAME(3:4),'(I2)') ISEC
          ELSE IF (ISEC.LE.999) THEN
            WRITE(CNAME(2:4),'(I3)') ISEC
          ELSE
            WRITE(CNAME,'(I4)') ISEC
          ENDIF
          CTYPE='SD'
          WRITE(99,'(1X,A4,1X,A2,2X,E12.6)') CNAME,CTYPE,X2-X1
          WRITE(98,'(4E16.6)')X1,X2,X2-X1,ZERO

          CALL SECTMAG_MAT(0.0D0,X2-X1,HMAT,VMAT)
          CALL UTIL_MATRIX_MULTIPLICATION(2,2,2,HMATT,HMAT,HMATT,WORK22)
          CALL UTIL_MATRIX_MULTIPLICATION(2,2,2,VMATT,VMAT,VMATT,WORK22)

        ENDIF !SPHI2.NE.SPHI1

        X1=X2
        Z1=X1

        PHI1=PHI2
        SPHI1=SPHI2

        IZERO=0

      IF (IX.LT.NCO) GOTO 1

      CLOSE(99)
      CLOSE(98)

      IF (HMATT(2,1).NE.0.0D0) THEN
        HFOC=-1.0D0/HMATT(2,1)
      ELSE
        HFOC=0.0D0
      ENDIF

      IF (VMATT(2,1).NE.0.0D0) THEN
        VFOC=-1.0D0/VMATT(2,1)
      ELSE
        VFOC=0.0D0
      ENDIF

      OPEN(UNIT=99,FILE='wave.smag',status='old')
      OPEN(UNIT=97,FILE='wave.strmag',status='new')

      WRITE(97,*)ICODE
      WRITE(97,*)CODE

      ISEC=0
      JSEC=0
      NSEC=0

      READ(99,'(1X,A4)',END=99)CNAME
      READ(99,'(1X,A4)',END=99)CNAME

11    READ(99,'(1X,A4)',END=99)CNAME

      NSEC=NSEC+1
      ISEC=ISEC+1

      CNAME10(ISEC)=CNAME

      IF (ISEC.EQ.10) THEN
        WRITE(97,'(10(1X,A4))')CNAME10
        CNAME10='NN  '
        ISEC=0
        JSEC=JSEC+1
      ENDIF

      GOTO 11

99    CONTINUE

      IF (ISEC.ne.0) THEN
        WRITE(97,'(10(1X,A4))')CNAME10
        JSEC=JSEC+1
      ENDIF

      CLOSE(97)
      CLOSE(99)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      SR SECTMAGS:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  '      Number of elements written to wave.smag and wave.xmag:',NSEC
      WRITE(LUNGFO,*)
     &'      Number of elements lines written to wave.strmag:',JSEC
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      Linear tranfer matrices:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      ',SNGL(HMATT(1,1)),SNGL(HMATT(1,2))
      WRITE(LUNGFO,*)'      ',SNGL(HMATT(2,1)),SNGL(HMATT(2,2))
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      ',SNGL(VMATT(1,1)),SNGL(VMATT(1,2))
      WRITE(LUNGFO,*)'      ',SNGL(VMATT(2,1)),SNGL(VMATT(2,2))
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      Focal lengths:'
      WRITE(LUNGFO,*)'      ',SNGL(HFOC),SNGL(VFOC)
      WRITE(LUNGFO,*)

      RETURN
      END
