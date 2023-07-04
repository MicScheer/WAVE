*CMZ :  4.00/13 12/11/2021  11.33.50  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/08 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.13/03 12/01/2000  16.31.33  by  Michael Scheer
*CMZ : 00.01/09 01/09/95  14.02.50  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  18.45.21  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.55.14  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.28  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE STOKSUMF
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

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.


      INTEGER IFREQ,ISTOK,IOBSV,IS,ICAL

      REAL*4 RAT(4)
      DOUBLE PRECISION DZ,DY
      DOUBLE PRECISION SUM(4),POL

      DATA ICAL/0/

      IF (IPIN.EQ.0) RETURN

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR STOKSUMF:'
      WRITE(LUNGFO,*)

         IF (IPINCIRC.NE.0.AND.ICAL.NE.1) THEN

             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** WARNING SR STOKSUMF ***'
             WRITE(LUNGFO,*)
     &'USE OF FLAG IPINCIRC NOT RECOMMENDED FOR THIS ROUTINE SINCE'
             WRITE(LUNGFO,*)
     &'SHAPE OF CIRCULARE PINHOLE IS TAKEN INTO ACCOUNT ONLY VERY ROUGHLY'
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
C            WRITE(6,*)
C            WRITE(6,*)
C            WRITE(6,*)'*** WARNING SR STOKSUMF ***'
C            WRITE(6,*)
C     &'USE OF FLAG IPINCIRC NOT RECOMMENDED FOR THIS ROUTINE SINCE'
C            WRITE(6,*)
C     &'SHAPE OF CIRCULARE PINHOLE IS TAKEN INTO ACCOUNT ONLY VERY ROUGHLY'
C            WRITE(6,*)
C            WRITE(6,*)
             ICAL=1
         ENDIF !IPINCIRC

      WRITE(LUNGFO,*)
     & '     Photon energy, S0,S1/S0,S2/S0,S3/S0 and polarization simply summed up'
      WRITE(LUNGFO,*)
     & '     (with emittance effects):'
      IF (IW_BLENF.NE.0) THEN
        WRITE(LUNGFO,*)'     (and ratio of spline and summation results)'
        WRITE(6,*)'     Ratio of spline and summation results:'
        write(6,*)
      ENDIF
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

         DO IFREQ=1,NFREQ

         DO ISTOK=1,4

             SUM(ISTOK)=0.0

         DO IOBSV=1,NOBSV

             IF (IPINCIRC.EQ.0) THEN

            IF (
     &           DABS(OBSV(2,IOBSV)-PINCEN(2))-PINH/2.D0.LT.1.D-10
     &                      .AND.
     &           DABS(OBSV(3,IOBSV)-PINCEN(3))-PINW/2.D0.LT.1.D-10
     &                      ) THEN

                  IF(DABS(
     &                           DABS(OBSV(3,IOBSV)-PINCEN(3))
     &                           -PINW/2.D0).LT.1.D-10) THEN

                DZ=OBSVDZ/2.D0

                  ELSE

                     DZ=OBSVDZ

                  ENDIF

                  IF(DABS(
     &                           DABS(OBSV(2,IOBSV)-PINCEN(2))
     &                           -PINH/2.D0).LT.1.D-10) THEN

                DY=OBSVDY/2.D0

                  ELSE

                    DY=OBSVDY

                  ENDIF

                  SUM(ISTOK)=SUM(ISTOK)+STOKESF(ISTOK,IOBSV+NOBSV*(IFREQ-1))*DZ*DY

            ENDIF   !OBSV

             ELSE    !IPINCIR

            IF (
     &                      (OBSV(2,IOBSV)-PINCEN(2))**2
     &                     +(OBSV(3,IOBSV)-PINCEN(3))**2
     &                       -PINR**2.LT.1.D-10) THEN

                  DZ=OBSVDZ
                  DY=OBSVDY
                  SUM(ISTOK)=SUM(ISTOK)+STOKESF(ISTOK,IOBSV+NOBSV*(IFREQ-1))*DZ*DY

                 ENDIF  !OBSV

             ENDIF   !IPINCIR

         ENDDO !NOBSV
         ENDDO !ISTOK


         IF (SUM(1).EQ.0.0) SUM(1)=-9999.
         POL=SQRT(SUM(2)**2+SUM(3)**2+
     &               SUM(4)**2)/SUM(1)

         IF (IUNIT.EQ.0)       !260194
     &      WRITE(LUNGFO,2584)SNGL(FREQ(IFREQ))
     &                ,SUM(1),(SUM(IS)/SUM(1),IS=2,4),POL
         IF (IUNIT.NE.0) !260194
     &      WRITE(LUNGFO,2584)SNGL(WELLEN(IFREQ))
     &                ,SUM(1),(SUM(IS)/SUM(1),IS=2,4),POL
2584     FORMAT('     ',6(1PE12.4))

        IF (IW_BLENF.NE.0) THEN
         DO IS=1,4
         RAT(IS)=9999.
         IF (SUM(IS).NE.0.D0) RAT(IS)=WSTOKESF(IS,IFREQ)/SUM(IS)
         ENDDO !IS
         WRITE(LUNGFO,2584)SNGL(FREQ(IFREQ)),(RAT(IS),IS=1,4)
         WRITE(6,2584)SNGL(FREQ(IFREQ)),(RAT(IS),IS=1,4)
        ENDIF   !IW_BLENF

         ENDDO !IFREQ

      RETURN
      END
