*CMZ :  2.61/02 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.53/05 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.52/14 22/12/2004  15.44.27  by  Michael Scheer
*CMZ :  2.52/13 16/12/2004  21.16.53  by  Michael Scheer
*-- Author :    Michael Scheer   15/12/2004
      SUBROUTINE BSECTMAGS(X,Y,Z,BX,BY,BZ,AX,AY,AZ)
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
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER NDIPP,ICAL
      PARAMETER(NDIPP=999)

      DOUBLE PRECISION DX,RHO,DPHI,RHODIP,
     &  PHIDIP(NDIPP),BYDIP(NDIPP),XIDIP(NDIPP),XEDIP(NDIPP),
     &  X,Y,Z,BX,BY,BZ,AX,AY,AZ

      INTEGER NSEC,NDIP,JCODE,LCODE,IDIPO,IDIP

      CHARACTER(4) CNAME
      CHARACTER(2) CTYPE
      CHARACTER(128) CODEJ,CODEL

      DATA ICAL/0/
      DATA IDIPO/1/

      IF (ICAL.EQ.0) THEN

        Y=Y
        Z=Z

        OPEN(UNIT=98,FILE='wave.xmag',STATUS='OLD')
        OPEN(UNIT=99,FILE='wave.smag',STATUS='OLD')

        READ(98,*)LCODE
        READ(98,'(A)')CODEL
        READ(99,*)JCODE
        READ(99,'(A)')CODEJ

        IF (LCODE.NE.JCODE) THEN
          WRITE(LUNGFO,*)
     &    '*** Error in BSECTMAGS: Files wave.smag and wave.xmag do not match'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** Programm WAVE aborted ***'
          WRITE(LUNGFO,*)
          WRITE(6,*)
     &    '*** Error in BSECTMAGS: Files wave.smag and wave.xmag do not match'
          WRITE(6,*)
          WRITE(6,*)'*** Programm WAVE aborted ***'
          WRITE(6,*)
          STOP
        ENDIF

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      Subroutine BSECTMAGS:'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      run number on file wave.smag:',JCODE
        WRITE(LUNGFO,*)'      comment on file wave.smag:'
        WRITE(LUNGFO,*)'      ',CODEJ
        WRITE(LUNGFO,*)

        NDIP=0

1       READ(99,'(1X,A4,1X,A2,2X,E12.6,2X)',END=9)
     &      CNAME,CTYPE,DX
        READ(98,*) DX

        NSEC=NSEC+1

        IF (NSEC.GT.999) THEN
          WRITE(LUNGFO,*)'*** Error in BSECTMAGS: Too many magnets'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** Programm WAVE aborted ***'
          WRITE(6,*)'*** Error in BSECTMAGS: Too many magnets'
          WRITE(6,*)
          WRITE(6,*)'*** Programm WAVE aborted ***'
          STOP
        ENDIF

        IF (CTYPE.EQ.'DI') THEN

          NDIP=NDIP+1

          BACKSPACE(99)
          BACKSPACE(98)

          READ(99,'(1X,A4,1X,A2,2X,E12.6,2X,E12.6)')
     &      CNAME,CTYPE,PHIDIP(NDIP),RHODIP
          READ(98,*)XIDIP(NDIP),XEDIP(NDIP),RHO,DPHI

          BYDIP(NDIP)=
     &      EMASSE1*DSQRT( (DMYGAMMA+1.0D0)*(DMYGAMMA-1.0D0))/CLIGHT1/RHO

          IF (DPHI.NE.PHIDIP(NDIP).OR.RHO.NE.RHODIP) THEN
            WRITE(LUNGFO,*)
     &   '*** Error in BSECTMAGS: Files wave.smag and wave.xmag do not match'
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** Programm WAVE aborted ***'
            WRITE(LUNGFO,*)
            WRITE(6,*)
     &   '*** Error in BSECTMAGS: Files wave.smag and wave.xmag do not match'
            WRITE(6,*)
            WRITE(6,*)'*** Programm WAVE aborted ***'
            WRITE(6,*)
            STOP
          ENDIF

        ENDIF !DI

        GOTO 1

9       CLOSE(99)
        CLOSE(98)

        IF (XSTART.EQ.9999.0D0) XSTART=XIDIP(1)
        IF (XSTOP.EQ.9999.0D0) XSTOP=XEDIP(NDIP)

        ICAL=1

      ENDIF !ICAL

      IF (X.LE.XEDIP(1)) THEN
        IDIPO=1
      ENDIF

      BY=0.0D0

      DO IDIP=IDIPO,NDIP
        IF (XIDIP(IDIP).LE.X.AND.XEDIP(IDIP).GE.X) THEN
          BY=BYDIP(IDIP)
          IDIPO=IDIP
          GOTO 99
        ENDIF
      ENDDO

99    BX=0.0D0
      BZ=0.0D0

      AX=0.0D0
      AY=0.0D0
      AZ=0.0D0

      RETURN
      END
