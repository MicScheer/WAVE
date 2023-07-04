*CMZ :  4.00/04 06/08/2019  10.50.42  by  Michael Scheer
*CMZ :  3.04/01 27/03/2018  12.56.47  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.61/01 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.52/16 17/01/2005  13.26.00  by  Michael Scheer
*CMZ :  2.48/04 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.48/03 03/03/2004  12.49.39  by  Michael Scheer
*CMZ :  2.46/02 07/03/2003  11.14.26  by  Michael Scheer
*CMZ :  2.41/13 22/08/2002  17.22.10  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.02  by  Michael Scheer
*CMZ :  2.40/03 14/03/2002  16.42.48  by  Michael Scheer
*CMZ :  2.40/02 14/03/2002  15.50.44  by  Michael Scheer
*CMZ :  2.40/00 11/03/2002  19.04.35  by  Michael Scheer
*CMZ :  2.39/00 14/12/2001  18.49.49  by  Michael Scheer
*CMZ :  2.37/00 09/11/2001  11.03.31  by  Michael Scheer
*CMZ :  2.20/01 17/01/2001  11.56.16  by  Michael Scheer
*CMZ :  2.16/08 31/10/2000  14.25.16  by  Michael Scheer
*CMZ :  2.16/05 24/08/2000  12.55.21  by  Michael Scheer
*CMZ :  2.16/04 20/07/2000  15.49.30  by  Michael Scheer
*CMZ :  1.00/00 30/09/97  11.30.53  by  Michael Scheer
*CMZ : 00.01/07 28/02/95  12.32.01  by  Michael Scheer
*CMZ : 00.01/05 01/02/95  15.24.06  by  Michael Scheer
*CMZ : 00.01/04 26/01/95  15.55.24  by  Michael Scheer
*CMZ : 00.00/07 24/05/94  09.48.34  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE USERASCII
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

C WRITE USER ASCII-FILES

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEND.

      INTEGER I,IFREQ,IY,IZ,IO,I1
      INTEGER IGETLASTCHAR,IGETFIRSTCHAR
      REAL A,B,C,D,E,F
      CHARACTER C1
      CHARACTER(16) C16
      CHARACTER(72) CLINE,FILE

      EXTERNAL IGETLASTCHAR,IGETFIRSTCHAR

C TRAJECTORY
      OPEN(UNIT=99,FILE='trajectory.wva')
      IF (IHISASCII.GT.0) THEN
        cline=""
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' Trajectory and mag. field X,Y,Z,BX,BY,BZ [m,T]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)') nco
        write(99,'(a)')cline
      ENDIF
      DO I=1,NCO
        A=WSXYZ(1,I)
        B=WSXYZ(2,I)
        C=WSXYZ(3,I)
        D=WBXYZ(1,I)
        E=WBXYZ(2,I)
        F=WBXYZ(3,I)
        WRITE(99,'(6(1PE14.6))') A,B,C,D,E,F
      ENDDO !NCO
      CLOSE(99)

C STOKES SINGLE POINT
      OPEN(UNIT=99,FILE='stokes_point.wva')
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' S0,S1,S2,S3 [ev, phot/s/BW/mm**2]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)')nfreq
        write(99,'(a)')cline
      ENDIF
      DO I=1,NFREQ
        A=FREQ(I)
        IF (ISTOKES.NE.0) THEN
          B=STOKEC(1,I)*1.E-6
          C=STOKEC(2,I)*1.E-6
          D=STOKEC(3,I)*1.E-6
          E=STOKEC(4,I)*1.E-6
        ELSE
          B=0.0
          C=0.0
          D=0.0
          E=0.0
        ENDIF
        WRITE(99,'(5(1PE14.6))') A,B,C,D,E
      ENDDO !NFREQ
      CLOSE(99)

C STOKES SINGLE POINT (FOLDED)
      OPEN(UNIT=99,FILE='stokes_point_emittance.wva')
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' S0,S1,S2,S3 (with emittance) [ev, phot/s/BW/mm**2]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)')nfreq
        write(99,'(a)')cline
      ENDIF
      DO I=1,NFREQ
        A=FREQ(I)
        IF (ISTOKES.NE.0.AND.IFOLD.NE.0) THEN
          B=STOKECF(1,I)*1.E-6
          C=STOKECF(2,I)*1.E-6
          D=STOKECF(3,I)*1.E-6
          E=STOKECF(4,I)*1.E-6
        ELSE
          B=0.0
          C=0.0
          D=0.0
          E=0.0
        ENDIF
        WRITE(99,'(5(1PE14.6))') A,B,C,D,E
      ENDDO !NFREQ
      CLOSE(99)

C STOKES SINGLE POINT (E-FOLDED)
      OPEN(UNIT=99,FILE='stokes_point_espread.wva')
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' S0,S1,S2,S3 (with beam energy spread) [ev, phot/s/BW/mm**2]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)')nfreq
        write(99,'(a)')cline
      ENDIF
      DO I=1,NFREQ
        A=FREQ(I)
        IF (ISTOKES.NE.0.AND.IEFOLD.NE.0) THEN
          B=STOKECE(1,I)*1.E-6
          C=STOKECE(2,I)*1.E-6
          D=STOKECE(3,I)*1.E-6
          E=STOKECE(4,I)*1.E-6
        ELSE
          B=0.0
          C=0.0
          D=0.0
          E=0.0
        ENDIF
        WRITE(99,'(5(1PE14.6))') A,B,C,D,E
      ENDDO !NFREQ
      CLOSE(99)

C STOKES SINGLE POINT (EF)
      OPEN(UNIT=99
     &  ,FILE='stokes_point_emittance_espread.wva')
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' S0,S1,S2,S3 (with emittance and energy spread) [ev, phot/s/BW/mm**2]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)')nfreq
        write(99,'(a)')cline
      ENDIF
      DO I=1,NFREQ
        A=FREQ(I)
        IF (ISTOKES.NE.0.AND.IFOLD.NE.0.AND.IEFOLD.NE.0) THEN
          B=STOKECEF(1,I)*1.E-6
          C=STOKECEF(2,I)*1.E-6
          D=STOKECEF(3,I)*1.E-6
          E=STOKECEF(4,I)*1.E-6
        ELSE
          B=0.0
          C=0.0
          D=0.0
          E=0.0
        ENDIF
        WRITE(99,'(5(1PE14.6))') A,B,C,D,E
      ENDDO !NFREQ
      CLOSE(99)

C BRILLIANCE
      OPEN(UNIT=99,FILE='stokes_brilliance.wva')
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' Brilliance of S0,S1,S2,S3 [ev, phot/s/BW/mm**2/mrad**2]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)')nfreq
        write(99,'(a)')cline
      ENDIF
      DO I=1,NFREQ
        A=FREQ(I)
        IF (IBRILL.NE.0) THEN
          B=BRILLC(1,I)*1.E-12
          C=BRILLC(2,I)*1.E-12
          D=BRILLC(3,I)*1.E-12
          E=BRILLC(4,I)*1.E-12
        ELSE
          B=0.0
          C=0.0
          D=0.0
          E=0.0
        ENDIF
        WRITE(99,'(5(1PE14.6))') A,B,C,D,E
      ENDDO !NFREQ
      CLOSE(99)

C BRILLIANCE (folded)
      OPEN(UNIT=99,FILE='stokes_brilliance_emittance.wva')
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' Brilliance of S0,S1,S2,S3 (with emittance) [ev, phot/s/BW/mm**2/mrad**2]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)')nfreq
        write(99,'(a)')cline
      ENDIF
      DO I=1,NFREQ
        A=FREQ(I)
        IF (IBRILL.NE.0.AND.IFOLD.NE.0) THEN
          B=BRILLCF(1,I)*1.E-12
          C=BRILLCF(2,I)*1.E-12
          D=BRILLCF(3,I)*1.E-12
          E=BRILLCF(4,I)*1.E-12
        ELSE
          B=0.0
          C=0.0
          D=0.0
          E=0.0
        ENDIF
        WRITE(99,'(5(1PE14.6))') A,B,C,D,E
      ENDDO !NFREQ
      CLOSE(99)

C BRILLIANCE (e-spread)
      OPEN(UNIT=99,FILE='stokes_brilliance_espread.wva')
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' Brilliance of S0,S1,S2,S3 (with energy spread) [ev, phot/s/BW/mm**2/mrad**2]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)')nfreq
        write(99,'(a)')cline
      ENDIF
      DO I=1,NFREQ
        A=FREQ(I)
        IF (IBRILL.NE.0.AND.IEFOLD.NE.0) THEN
          B=BRILLCE(1,I)*1.E-12
          C=BRILLCE(2,I)*1.E-12
          D=BRILLCE(3,I)*1.E-12
          E=BRILLCE(4,I)*1.E-12
        ELSE
          B=0.0
          C=0.0
          D=0.0
          E=0.0
        ENDIF
        WRITE(99,'(5(1PE14.6))') A,B,C,D,E
      ENDDO !NFREQ
      CLOSE(99)

C BRILLIANCE (emittance, e-spread)
      OPEN(UNIT=99,FILE='stokes_brilliance_emittance_espread.wva')

      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' Brilliance of S0,S1,S2,S3 (w. emit. and e-spread) [ev, phot/s/BW/mm**2/mrad**2]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)')nfreq
        write(99,'(a)')cline
      ENDIF
      DO I=1,NFREQ
        A=FREQ(I)
        IF (IBRILL.NE.0.AND.IEFOLD.NE.0) THEN
          B=BRILLCEF(1,I)*1.E-12
          C=BRILLCEF(2,I)*1.E-12
          D=BRILLCEF(3,I)*1.E-12
          E=BRILLCEF(4,I)*1.E-12
        ELSE
          B=0.0
          C=0.0
          D=0.0
          E=0.0
        ENDIF
        WRITE(99,'(5(1PE14.6))') A,B,C,D,E
      ENDDO !NFREQ
      CLOSE(99)

C STOKES PINHOLE
      OPEN(UNIT=99,FILE='stokes_pinhole.wva')

      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' S0,S1,S2,S3 [ev, phot/s/BW]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)')nfreq
        write(99,'(a)')cline
      ENDIF
      DO I=1,NFREQ
        A=FREQ(I)
        IF (ISTOKES.NE.0) THEN
          B=WSTOKES(1,I)
          C=WSTOKES(2,I)
          D=WSTOKES(3,I)
          E=WSTOKES(4,I)
        ELSE
          B=0.0
          C=0.0
          D=0.0
          E=0.0
        ENDIF
        WRITE(99,'(5(1PE14.6))') A,B,C,D,E
      ENDDO !NFREQ
      CLOSE(99)

C STOKES (FOLDED)
      OPEN(UNIT=99,FILE='stokes_pinhole_emittance.wva')
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' S0,S1,S2,S3 (with emittance) [ev, phot/s/BW]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)')nfreq
        write(99,'(a)')cline
      ENDIF
      DO I=1,NFREQ
        A=FREQ(I)
        IF (ISTOKES.NE.0.AND.IFOLD.NE.0) THEN
          B=WSTOKESF(1,I)
          C=WSTOKESF(2,I)
          D=WSTOKESF(3,I)
          E=WSTOKESF(4,I)
        ELSE
          B=0.0
          C=0.0
          D=0.0
          E=0.0
        ENDIF
        WRITE(99,'(5(1PE14.6))') A,B,C,D,E
      ENDDO !NFREQ
      CLOSE(99)

C STOKES (E-FOLDED)
      OPEN(UNIT=99,FILE='stokes_pinhole_espread.wva')
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' S0,S1,S2,S3 (with beam energy spread) [ev, phot/s/BW]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)')nfreq
        write(99,'(a)')cline
      ENDIF
      DO I=1,NFREQ
        A=FREQ(I)
        IF (ISTOKES.NE.0.AND.IEFOLD.NE.0) THEN
          B=WSTOKESE(1,I)
          C=WSTOKESE(2,I)
          D=WSTOKESE(3,I)
          E=WSTOKESE(4,I)
        ELSE
          B=0.0
          C=0.0
          D=0.0
          E=0.0
        ENDIF
        WRITE(99,'(5(1PE14.6))') A,B,C,D,E
      ENDDO !NFREQ
      CLOSE(99)

C STOKES (EF)
      OPEN(UNIT=99
     &  ,FILE='stokes_pinhole_emittance_espread.wva')
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' S0,S1,S2,S3 (with emittance and energy spread) [ev, phot/s/BW]'
        cline(3:72)=""
        WRITE(CLINE(3:10),'(I8)')nfreq
        write(99,'(a)')cline
      ENDIF
      DO I=1,NFREQ
        A=FREQ(I)
        IF (ISTOKES.NE.0.AND.IFOLD.NE.0.AND.IEFOLD.NE.0) THEN
          B=WSTOKESEF(1,I)
          C=WSTOKESEF(2,I)
          D=WSTOKESEF(3,I)
          E=WSTOKESEF(4,I)
        ELSE
          B=0.0
          C=0.0
          D=0.0
          E=0.0
        ENDIF
        WRITE(99,'(5(1PE14.6))') A,B,C,D,E
      ENDDO !NFREQ
      CLOSE(99)

      if (ipin.ne.0) then

C STOKES DISTRIBUTION IN PINHOLE
      DO IFREQ=1,NFREQ
        WRITE(C16,*)IFREQ
        I=IGETLASTCHAR(1,16,C16,C1)
        I1=IGETFIRSTCHAR(1,16,C16,C1)
        IF (IFREQ.EQ.1) THEN
          FILE='stokes_dist_'//C16(I1:I)//'st_energy.wva'
        ELSE IF (IFREQ.EQ.2) THEN
          FILE='stokes_dist_'//C16(I1:I)//'nd_energy.wva'
        ELSE IF (IFREQ.EQ.3) THEN
          FILE='stokes_dist_'//C16(I1:I)//'rd_energy.wva'
        ELSE
          FILE='stokes_dist_'//C16(I1:I)//'th_energy.wva'
        ENDIF
        OPEN(UNIT=99,FILE=FILE)
        IF (IHISASCII.GT.0) THEN
          CLINE(1:1)=CHISASCII
          CLINE(2:2)=' '
          WRITE(CLINE(3:10),'(I8)')ICODE
          CLINE=CLINE(1:10)//' | '//CODE
          write(99,'(a)')cline
          write(99,'(a)')chisascii //
     &      ' Z,Y,S0,S1,S2,S3 [mm, phot/s/BW/mm**2]'
          cline(3:72)=""
          write(cline(3:72),*)chisascii //' photon energy: ',SNGL(FREQ(IFREQ))
          write(99,'(a)')cline
        ENDIF

        DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY
          DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ

            IO=(IY-1)*NOBSVZ+IZ
            I=IO+NOBSV*(IFREQ-1)
            A=OBSV(3,IO)*1.E3
            B=OBSV(2,IO)*1.E3
            IF (ISTOKES.NE.0) THEN
              C=STOKES(1,I)*1.E-6
              D=STOKES(2,I)*1.E-6
              E=STOKES(3,I)*1.E-6
              F=STOKES(4,I)*1.E-6
            ELSE
              C=0.0
              D=0.0
              E=0.0
              F=0.0
            ENDIF
            WRITE(99,'(6(1PE14.6))') A,B,C,D,E,F
          ENDDO   !IZ
        ENDDO  !IY
        CLOSE(99)
      ENDDO !NFREQ

C STOKES DISTRIBUTION IN PINHOLE (FOLDED)
      DO IFREQ=1,NFREQ
        WRITE(C16,*)IFREQ
        I=IGETLASTCHAR(1,16,C16,C1)
        I1=IGETFIRSTCHAR(1,16,C16,C1)
        IF (IFREQ.EQ.1) THEN
          FILE='stokes_dist_emittance_'//C16(I1:I)//'st_energy.wva'
        ELSE IF (IFREQ.EQ.2) THEN
          FILE='stokes_dist_emittance_'//C16(I1:I)//'nd_energy.wva'
        ELSE IF (IFREQ.EQ.3) THEN
          FILE='stokes_dist_emittance_'//C16(I1:I)//'rd_energy.wva'
        ELSE
          FILE='stokes_dist_emittance_'//C16(I1:I)//'th_energy.wva'
        ENDIF
        OPEN(UNIT=99,FILE=FILE)
        IF (IHISASCII.GT.0) THEN
          CLINE(1:1)=CHISASCII
          CLINE(2:2)=' '
          WRITE(CLINE(3:10),'(I8)')ICODE
          CLINE=CLINE(1:10)//' | '//CODE
          write(99,'(a)')cline
          write(99,'(a)')chisascii //
     &      ' Z,Y,S0,S1,S2,S3 (with emittance) [mm, phot/s/BW/mm**2]'
          cline(3:72)=""
          write(cline(3:72),*)chisascii //' photon energy: ',SNGL(FREQ(IFREQ))
          write(99,'(a)')cline
        ENDIF

        DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY
          DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ
            IO=(IY-1)*NOBSVZ+IZ
            I=IO+NOBSV*(IFREQ-1)
            A=OBSV(3,IO)*1.E3
            B=OBSV(2,IO)*1.E3
            IF (ISTOKES.NE.0.AND.IFOLD.NE.0) THEN
              C=STOKESF(1,I)*1.E-6
              D=STOKESF(2,I)*1.E-6
              E=STOKESF(3,I)*1.E-6
              F=STOKESF(4,I)*1.E-6
            ELSE
              C=0.0
              D=0.0
              E=0.0
              F=0.0
            ENDIF
            WRITE(99,'(6(1PE14.6))') A,B,C,D,E,F
          ENDDO   !IZ
        ENDDO  !IY
        CLOSE(99)
      ENDDO !NFREQ

C IRRADIATED POWER DISTRIBUTION IN PINHOLE
      FILE='irradiated_power_dist.wva'
      OPEN(UNIT=99,FILE=FILE)
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' Z,Y,totally irradiated power-density [mm, W/mm**2]'
      ENDIF

      DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY
        DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ
          IO=(IY-1)*NOBSVZ+IZ
          A=OBSV(3,IO)*1.E3
          B=OBSV(2,IO)*1.E3
          C=SPECPOWT(IO)*1.E-6
          WRITE(99,'(3(1PE14.6))') A,B,C
        ENDDO  !IZ
      ENDDO !IY
      CLOSE(99)

C POWER DISTRIBUTION IN PINHOLE
      FILE='power_dist.wva'
      OPEN(UNIT=99,FILE=FILE)
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &' Z,Y,power-density (from spectrum) [mm, W/mm**2]'
      ENDIF

      DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY
        DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ
          IO=(IY-1)*NOBSVZ+IZ
          A=OBSV(3,IO)*1.E3
          B=OBSV(2,IO)*1.E3
          C=SPECTOTI(IO)*1.E-6
          WRITE(99,'(3(1PE14.6))') A,B,C
        ENDDO  !IZ
      ENDDO !IY
      CLOSE(99)

C POWER DISTRIBUTION IN PINHOLE (folded)
      FILE='power_dist_emittance.wva'
      OPEN(UNIT=99,FILE=FILE)
      IF (IHISASCII.GT.0) THEN
        CLINE(1:1)=CHISASCII
        CLINE(2:2)=' '
        WRITE(CLINE(3:10),'(I8)')ICODE
        CLINE=CLINE(1:10)//' | '//CODE
        write(99,'(a)')cline
        write(99,'(a)')chisascii //
     &    ' Z,Y,power-density (from spectrum, with emittance) [mm, W/mm**2]'
      ENDIF

      DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY
        DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ
          IO=(IY-1)*NOBSVZ+IZ
          A=OBSV(3,IO)*1.E3
          B=OBSV(2,IO)*1.E3
          IF (IFOLD.NE.0) THEN
            C=SPECTOTIF(IO)*1.E-6
          ELSE
            C=0.0
          ENDIF
          WRITE(99,'(3(1PE14.6))') A,B,C
        ENDDO  !IZ
      ENDDO !IY
      CLOSE(99)

C STOKES DISTRIBUTION IN PINHOLE
      DO IFREQ=1,NFREQ

        WRITE(C16,*)IFREQ
        I=IGETLASTCHAR(1,16,C16,C1)
        I1=IGETFIRSTCHAR(1,16,C16,C1)
        IF (IFREQ.EQ.1) THEN
          FILE='stokes_dist_espread_'//C16(I1:I)//'st_energy.wva'
        ELSE IF (IFREQ.EQ.2) THEN
          FILE='stokes_dist_espread_'//C16(I1:I)//'nd_energy.wva'
        ELSE IF (IFREQ.EQ.3) THEN
          FILE='stokes_dist_espread_'//C16(I1:I)//'rd_energy.wva'
        ELSE
          FILE='stokes_dist_espread_'//C16(I1:I)//'th_energy.wva'
        ENDIF
        OPEN(UNIT=99,FILE=FILE)
        IF (IEFOLD.NE.0) THEN
          IF (IHISASCII.GT.0) THEN
            CLINE(1:1)=CHISASCII
            CLINE(2:2)=' '
            WRITE(CLINE(3:10),'(I8)')ICODE
            CLINE=CLINE(1:10)//' | '//CODE
            write(99,'(a)')cline
            write(99,'(a)')chisascii //
     &        ' Z,Y,S0,S1,S2,S3 (with e-spread) [mm, phot/s/BW/mm**2]'
            cline(3:72)=""
            write(cline(3:72),*)chisascii //' photon energy: ',SNGL(FREQ(IFREQ))
            write(99,'(a)')cline
          ENDIF

          DO IO=1,NOBSV

            I=IO+NOBSV*(IFREQ-1)

            A=OBSV(3,IO)*1.E3
            B=OBSV(2,IO)*1.E3

            C=STOKESE(1,I)*1.E-6
            D=STOKESE(2,I)*1.E-6
            E=STOKESE(3,I)*1.E-6
            F=STOKESE(4,I)*1.E-6

            WRITE(99,'(6(1PE14.6))') A,B,C,D,E,F

          ENDDO   !NOBSV

        ENDIF !IEFOLD

        CLOSE(99)

      ENDDO !NFREQ

      DO IFREQ=1,NFREQ

        WRITE(C16,*)IFREQ
        I=IGETLASTCHAR(1,16,C16,C1)
        I1=IGETFIRSTCHAR(1,16,C16,C1)
        IF (IFREQ.EQ.1) THEN
          FILE='stokes_dist_emittance_espread_'//C16(I1:I)//'st_energy.wva'
        ELSE IF (IFREQ.EQ.2) THEN
          FILE='stokes_dist_emittance_espread_'//C16(I1:I)//'nd_energy.wva'
        ELSE IF (IFREQ.EQ.3) THEN
          FILE='stokes_dist_emittance_espread_'//C16(I1:I)//'rd_energy.wva'
        ELSE
          FILE='stokes_dist_emittance_espread_'//C16(I1:I)//'th_energy.wva'
        ENDIF
        OPEN(UNIT=99,FILE=FILE)

        IF (IFOLD.NE.0.AND.IEFOLD.NE.0) THEN

          IF (IHISASCII.GT.0) THEN
            CLINE(1:1)=CHISASCII
            CLINE(2:2)=' '
            WRITE(CLINE(3:10),'(I8)')ICODE
            CLINE=CLINE(1:10)//' | '//CODE
            write(99,'(a)')cline
            write(99,'(a)')chisascii //
     &        ' Z,Y,S0,S1,S2,S3 (with emittance and e-spread) [mm, phot/s/BW/mm**2]'
            cline(3:72)=""
            write(cline(3:72),*)chisascii //' photon energy: ',SNGL(FREQ(IFREQ))
            write(99,'(a)')cline
          ENDIF

          DO IO=1,NOBSV

            I=IO+NOBSV*(IFREQ-1)

            A=OBSV(3,IO)*1.E3
            B=OBSV(2,IO)*1.E3

            C=STOKESEF(1,I)*1.E-6
            D=STOKESEF(2,I)*1.E-6
            E=STOKESEF(3,I)*1.E-6
            F=STOKESEF(4,I)*1.E-6

            WRITE(99,'(6(1PE14.6))') A,B,C,D,E,F

          ENDDO   !NOBSV

        ENDIF !IFOLD, IEFOLD

        CLOSE(99)

      ENDDO !NFREQ

      endif !(ipin.ne.0 then

      RETURN
      END
