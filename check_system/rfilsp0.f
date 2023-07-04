*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.41/10 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.45  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.35  by  Michael Scheer
*CMZ :  2.13/03 10/01/2000  17.17.14  by  Michael Scheer
*CMZ :  2.00/00 05/01/99  16.19.21  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.47.16  by  Michael Scheer
*CMZ : 00.02/04 24/02/97  12.37.49  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.18.53  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.53.53  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.17  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE RFILSP0
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
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEND.

      CHARACTER(60) CODSP0
      INTEGER ICODSP0,IO,IX,IFR,IS,N2POWY,N2POWZ,IY,IZ,IOBSV,
     &          NOBSVZN,NOBSVYN,NOBSVN,IOO

         OPEN(UNIT=LUNSP0,FILE=FILESP0,STATUS='OLD')

         READ(LUNSP0,'(I12,A60)')ICODSP0,CODSP0
         READ(LUNSP0,*)

         IF (IUSEM .NE. 0) THEN


         DO IS=1,NSOURCE
         DO IO=1,NOBSV
         DO IFR=1,NFREQ
             SPEC(IS+NSOURCE*(IO-1+NOBSV*(IFR-1)))=0.0
         ENDDO
         ENDDO
         ENDDO

         READ(LUNSP0,*)NSOURCE,NOBSVN,NFREQ,IFREQ2P
         READ(LUNSP0,*)NOBSVZN,NOBSVYN,MOBSVZ,MOBSVY
         READ(LUNSP0,*)MEDGEZ,MEDGEY,MMEDGEZ,MMEDGEY
         READ(LUNSP0,*)
         READ(LUNSP0,*)PINW,PINH,PINR
         READ(LUNSP0,*)OBSVDZ,OBSVDY
         READ(LUNSP0,*)

           IF (NOBSVZN/2*2.NE.NOBSVZN) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** ERROR IN RFILSP0 ***'
             WRITE(LUNGFO,*)
     &'FLAG IUSEM IS SET, BUT NUMBER OF HORIZONTAL GRID POINTS IS NOT EVEN!'
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)'*** ERROR IN RFILSP0 ***'
             WRITE(6,*)
     &'FLAG IUSEM IS SET, BUT NUMBER OF HORIZONTAL GRID POINTS IS NOT EVEN!'
             WRITE(6,*)
             WRITE(6,*)
             STOP
           ENDIF  !NOBSVZN
           IF (NOBSVYN/2*2.NE.NOBSVYN) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** ERROR IN RFILSP0 ***'
             WRITE(LUNGFO,*)
     &'FLAG IUSEM IS SET, BUT NUMBER OF VERTICAL GRID POINTS IS NOT EVEN!'
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)'*** ERROR IN RFILSP0 ***'
             WRITE(6,*)
     &'FLAG IUSEM IS SET, BUT NUMBER OF VERTICAL GRID POINTS IS NOT EVEN!'
             WRITE(6,*)
             WRITE(6,*)
             STOP
           ENDIF  !NOBSVYN

           N2POWZ=NINT(ALOG(FLOAT(NOBSVZN-1))/ALOG(2.))
           IF(NOBSVZN .GT. 2**N2POWZ) N2POWZ=N2POWZ+1
           NOBSVZ=2**N2POWZ

           N2POWY=NINT(ALOG(FLOAT(NOBSVYN-1))/ALOG(2.))
           IF(NOBSVYN .GT. 2**N2POWY) N2POWY=N2POWY+1
           NOBSVY=2**N2POWY

           NOBSV=NOBSVZ*NOBSVY

           IF (NOBSVZ.GT.NDOBSVZP) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** ERROR IN RFILSP0 ***'
             WRITE(LUNGFO,*)
     &'DIMENSION EXCEEDED, INCREASE PARAMETER NDOBSVZP IN CMPARA.CMN'
             WRITE(LUNGFO,*)'MUST BE AT LEAST:',NOBSVZ
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)
             WRITE(6,*)'*** ERROR IN RFILSP0 ***'
             WRITE(6,*)
     &'DIMENSION EXCEEDED, INCREASE PARAMETER NDOBSVZP IN CMPARA.CMN'
             WRITE(6,*)'MUST BE AT LEAST:',NOBSVZ
             WRITE(6,*)
             WRITE(6,*)
             STOP
           ENDIF      !NOBVZP

           IF (NOBSVY.GT.NDOBSVYP) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** ERROR IN RFILSP0 ***'
             WRITE(LUNGFO,*)
     &'DIMENSION EXCEEDED, INCREASE PARAMETER NDOBSVYP IN CMPARA.CMN'
             WRITE(LUNGFO,*)'MUST BE AT LEAST:',NOBSVY
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)
             WRITE(6,*)'*** ERROR IN RFILSP0 ***'
             WRITE(6,*)
     &'DIMENSION EXCEEDED, INCREASE PARAMETER NDOBSVYP IN CMPARA.CMN'
             WRITE(6,*)'MUST BE AT LEAST:',NOBSVY
             WRITE(6,*)
             WRITE(6,*)
             STOP
           ENDIF      !NOBVZP

             MEDGEZ=(NOBSVZ-MOBSVZ)/2
             MEDGEY=(NOBSVY-MOBSVY)/2
             MMEDGEZ=0
             MMEDGEY=0

         READ(LUNSP0,*)(OBSVZ(IO),
     &              IO=(NOBSVZ-NOBSVZN)/2+1,(NOBSVZ-NOBSVZN)/2+NOBSVZN)
         READ(LUNSP0,*)
         READ(LUNSP0,*)(OBSVY(IO),
     &              IO=(NOBSVY-NOBSVYN)/2+1,(NOBSVY-NOBSVYN)/2+NOBSVYN)
         READ(LUNSP0,*)

         DO IO=1,NOBSVN

            IZ=MOD(IO-1,NOBSVZN)+1
            IY=(IO-1)/NOBSVZN+1
            IOO=((NOBSVY-NOBSVYN)/2+IY-1)*NOBSVZ+((NOBSVZ-NOBSVZN)/2+IZ)
            READ(LUNSP0,*)(OBSV(IX,IOO),IX=1,3)

         DO IFR=1,NFREQ

            READ(LUNSP0,*)FREQ(IFR)
     &      ,(SPEC(IS+NSOURCE*(IOO-1+NOBSV*(IFR-1))),IS=1,NSOURCE)
     &                        ,SPECTOT(IOO+NOBSV*(IFR-1))

         ENDDO
         ENDDO

         IOBSV=0
         DO IY=1,NOBSVY
         DO IZ=1,NOBSVZ
             IOBSV=IOBSV+1
             OBSVZ(IZ)=-DFLOAT(NOBSVZ-1)/2.*OBSVDZ+(IZ-1)*OBSVDZ
             OBSVY(IY)=-DFLOAT(NOBSVY-1)/2.*OBSVDY+(IY-1)*OBSVDY
             OBSV(3,IOBSV)=-DFLOAT(NOBSVZ-1)/2.*OBSVDZ+(IZ-1)*OBSVDZ
             OBSV(2,IOBSV)=-DFLOAT(NOBSVY-1)/2.*OBSVDY+(IY-1)*OBSVDY
         ENDDO
         ENDDO

         ELSE  !IUSEM

         READ(LUNSP0,*)NSOURCE,NOBSV,NFREQ,IFREQ2P
         READ(LUNSP0,*)NOBSVZ,NOBSVY,MOBSVZ,MOBSVY
         READ(LUNSP0,*)MEDGEZ,MEDGEY,MMEDGEZ,MMEDGEY
         READ(LUNSP0,*)
         READ(LUNSP0,*)PINW,PINH,PINR
         READ(LUNSP0,*)OBSVDZ,OBSVDY
         READ(LUNSP0,*)

         READ(LUNSP0,*)(OBSVZ(IO),IO=1,NOBSVZ)
         READ(LUNSP0,*)
         READ(LUNSP0,*)(OBSVY(IO),IO=1,NOBSVY)
         READ(LUNSP0,*)

         DO IO=1,NOBSV
            READ(LUNSP0,*)(OBSV(IX,IO),IX=1,3)
         DO IFR=1,NFREQ
            READ(LUNSP0,*)FREQ(IFR)
     &      ,(SPEC(IS+NSOURCE*(IO-1+NOBSV*(IFR-1))),IS=1,NSOURCE)
     &                        ,SPECTOT(IO+NOBSV*(IFR-1))

         ENDDO
         ENDDO

         ENDIF !IUSEM

         CLOSE(LUNSP0)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'*** MESSAGE SR RFILSP0 ***'
      WRITE(LUNGFO,*)'DATA OF SPECTRUM CALCULATIONS READ FROM FILE:'
      WRITE(LUNGFO,*)FILESP0
      WRITE(LUNGFO,*)
     &  'OLD DATA OF NSOURCE,NOBSV,NFREQ,OBSV,FREQ,SPEC,SPECTOT ETC. OVERWRITTEN'
      WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
      WRITE(LUNGFO,*)'JOBNUMBER, USER COMMENT:'
      WRITE(LUNGFO,*)ICODSP0,'; ',CODSP0
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

      IF (IUSEM.NE.0) THEN

C--- CHECK IF NUMBER OF GRID POINTS ARE POWERS OF 2

      N2POWZ=NINT(ALOG(FLOAT(NOBSVZ-1))/ALOG(2.))

      IF (NOBSVZ.NE.2**N2POWZ) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN RFILSP0 ***'
          WRITE(LUNGFO,*)
     &      'Number of horizontal observation points not power of 2'
          WRITE(LUNGFO,*)'Check input file WAVE.IN'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN RFILSP0 ***'
          WRITE(6,*)
     &      'Number of horizontal observation points not power of 2'
          WRITE(6,*)'Check input file WAVE.IN'
          WRITE(6,*)
          STOP
      ENDIF

      N2POWY=NINT(ALOG(FLOAT(NOBSVY-1))/ALOG(2.))

      IF (NOBSVY.NE.2**N2POWY) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN RFILSP0 ***'
          WRITE(LUNGFO,*)
     &      'Number of horizontal observation points not power of 2'
          WRITE(LUNGFO,*)'Check input file WAVE.IN'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN RFILSP0 ***'
          WRITE(6,*)
     &      'Number of horizontal observation points not power of 2'
          WRITE(6,*)'Check input file WAVE.IN'
          WRITE(6,*)
          STOP
      ENDIF
      ENDIF !IUSEM

      RETURN
      END
