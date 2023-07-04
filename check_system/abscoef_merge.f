*CMZ :  2.70/12 01/03/2013  15.45.11  by  Michael Scheer
*CMZ :  2.41/10 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.37/03 23/11/2001  18.02.57  by  Michael Scheer
*CMZ :  2.36/00 07/11/2001  14.24.11  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ : 00.01/04 28/11/94  18.52.29  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  13.57.04  by  Michael Scheer
*CMZ : 00.01/01 21/06/94  15.22.20  by  Michael Scheer
*-- Author :    Michael Scheer   20/06/94

      SUBROUTINE ABSCOEF_MERGE

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

C Merges files with different absorption coefficients
C The routine reads a list of file names given in the
C the file FILEAM and writes the results to FILEAMO.
C ABSTHI (namelist $FREQN) and FILESP0 (namelist $MYFILES) are
C superseded by THICKTOT and FILEAMO

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEND.

      CHARACTER(65) CDUM,CFILE

      INTEGER NFILES,IFILES,LUNDUM,IERR,IEDGE
      INTEGER IETOT,NETOT,METOT,NENE,IENE,IKNOWN

      INTEGER NEBUFFP
      PARAMETER(NEBUFFP=1000)
      DOUBLE PRECISION EBUFF(NEBUFFP),ABUFF(NEBUFFP),EPHOT
      DOUBLE PRECISION EBUFFT(NEBUFFP)
      DOUBLE PRECISION THICK,DENS,COEFTOT(NEBUFFP),THICKTOT,DENSTOT,DENTHIT
      DOUBLE PRECISION AMU

      DATA LUNDUM /10/

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR ABSCOEF_MERGE:'
      WRITE(LUNGFO,*)'     ================='
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     List of files to be merged:'
      WRITE(LUNGFO,*)

C--- Open list of files to be merged

      OPEN(UNIT=LUNAM,FILE=FILEAM,STATUS='OLD')

      NFILES=0
100   CONTINUE

      READ(LUNAM,*,END=99) THICK,CFILE
      NFILES=NFILES+1
      WRITE(LUNGFO,*)'     ',CFILE
      WRITE(LUNGFO,*)'     Thickness [m]:',THICK
      WRITE(LUNGFO,*)
      GOTO 100

99    CONTINUE
      REWIND(LUNAM)


      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

C--  Read files and merge photon energies

         READ(LUNAM,*) THICK,CFILE
         OPEN(UNIT=LUNDUM,FILE=CFILE,STATUS='OLD')
            READ(LUNDUM,'(A65)')CDUM
            READ(LUNDUM,*)DENS
            READ(LUNDUM,*)NENE
            DO IENE=1,NENE
                READ(LUNDUM,*)EBUFFT(IENE)
            ENDDO !NENE
         CLOSE(LUNDUM)
         NETOT=NENE
         METOT=NETOT

      DO IFILES=2,NFILES

         READ(LUNAM,*) THICK,CFILE
         OPEN(UNIT=LUNDUM,FILE=CFILE,STATUS='OLD')
            READ(LUNDUM,'(A65)')CDUM
            READ(LUNDUM,*)DENS
            READ(LUNDUM,*)NENE
          DO IENE=1,NENE

                READ(LUNDUM,*)EPHOT

C-   Check if energy is already in buffer
C    Be careful with absorption edges;
c    edges are defined by two identical photon energies.
c    Check carefully if files contain identical elements


               IKNOWN=0
               DO IETOT=1,NETOT
                 IF (EBUFFT(IETOT).EQ.EPHOT) THEN
                   IKNOWN=1
                   GOTO 200
                 ENDIF
               ENDDO !NETOT
200            CONTINUE

               IF (IKNOWN.EQ.0) THEN

                   METOT=METOT+1  !DO NOT UPDATE NETOT HERE, BECAUSE OF EDGES

                   IF (METOT.GT.NEBUFFP) THEN
                     WRITE(6,*)
                     WRITE(6,*)
     & '*** ERROR IN ABSCOEF_MERGE: DIMENSION EXCEEDED, INCREASE NEBUFFP IN THIS ROUTINE ***'
                     WRITE(6,*)
                     WRITE(LUNGFO,*)
                     WRITE(LUNGFO,*)
     & '*** ERROR IN ABSCOEF_MERGE: DIMENSION EXCEEDED, INCREASE NEBUFFP IN THIS ROUTINE ***'
                     WRITE(LUNGFO,*)
                     STOP
                   ENDIF   !METOT
                   EBUFFT(METOT)=EPHOT
                 ENDIF    !IKNOWN

               ENDDO   !IETOT

            NETOT=METOT

         CLOSE(LUNDUM)

      ENDDO !IFILES
      REWIND(LUNAM)

C-  Sort buffer EBUF

            CALL UTIL_SORT(NETOT,EBUFFT)

C--- Do the merging

      THICKTOT=0.0
      DENSTOT=0.0
      DENTHIT=0.0
      DO IETOT=1,NETOT
             COEFTOT(IETOT)=0.0
      ENDDO   !NETOT

      DO IFILES=1,NFILES

      READ(LUNAM,*) THICK,CFILE

      OPEN(UNIT=LUNDUM,FILE=CFILE,STATUS='OLD')
         READ(LUNDUM,'(A65)')CDUM
         READ(LUNDUM,*)DENS
         READ(LUNDUM,*)NENE
         THICKTOT=THICKTOT+THICK
         DENTHIT=DENTHIT+DENS*THICK
         DO IENE=1,NENE
             READ(LUNDUM,*)EBUFF(IENE),ABUFF(IENE)
         ENDDO !IENE
      CLOSE(LUNDUM)

C    Be careful with absorption edges;
c    edges are defined by two identical photon energies.
c    Check carefully if files contain identical elements

         DO IETOT=1,NETOT
             EPHOT=EBUFFT(IETOT)
           IF (IFILTER.EQ.-1) THEN
               CALL ABSNOSPLI(EBUFF,ABUFF,NENE,EPHOT,AMU,IERR,1)
           ELSE
               CALL ABSNOSPLI(EBUFF,ABUFF,NENE,EPHOT,AMU,IERR,-1)
           ENDIF
             IF (IERR.NE.0) THEN
               WRITE(6,*)
               WRITE(6,*)'*** ERROR IN ABSCOEF_MERGE ***'
               WRITE(6,*)'INTERPOLATION OF COEFFICIENT FAILED'
               WRITE(6,*)'FOR PHOTON ENERGY:',EPHOT
               WRITE(6,*)'AND FILE:',CFILE
               WRITE(6,*)'CHECK THIS FILE'
               WRITE(6,*)
               WRITE(LUNGFO,*)
               WRITE(LUNGFO,*)'*** ERROR IN ABSCOEF_MERGE ***'
               WRITE(LUNGFO,*)'INTERPOLATION OF COEFFICIENT FAILED'
               WRITE(LUNGFO,*)'FOR PHOTON ENERGY:',EPHOT
               WRITE(LUNGFO,*)'AND FILE:',CFILE
               WRITE(LUNGFO,*)'CHECK THIS FILE'
               WRITE(LUNGFO,*)
               STOP
             ENDIF
C- Edge
            IEDGE=0

            IF (IETOT.LT.NETOT) THEN
               IF (EBUFFT(IETOT).EQ.EBUFFT(IETOT+1)) THEN
               DO IENE=1,NENE
            IF (EBUFF(IENE).EQ.EPHOT) IEDGE=IEDGE+1
               ENDDO
               IF (IEDGE.EQ.2) THEN
                       IEDGE=0
                  DO IENE=1,NENE
               IF (EBUFF(IENE).EQ.EPHOT) THEN
                   IEDGE=IENE
                   GOTO 17
               ENDIF
                  ENDDO
17                     AMU=ABUFF(IENE)
                    ENDIF !(IEDGE.EQ.2)
               ENDIF   !EDGE
               ENDIF   !IETOT.LT.NETOT

            IF (IETOT.GT.1) THEN
               IF (EBUFFT(IETOT).EQ.EBUFFT(IETOT-1)) THEN
               DO IENE=1,NENE
            IF (EBUFF(IENE).EQ.EPHOT) IEDGE=IEDGE+1
               ENDDO
               IF (IEDGE.EQ.2) THEN
                       IEDGE=0
                  DO IENE=1,NENE
               IF (EBUFF(IENE).EQ.EPHOT) THEN
                   IEDGE=IENE
                   GOTO 18
               ENDIF
                  ENDDO
18                AMU=ABUFF(IENE+1)
                    ENDIF !(IEDGE.EQ.2)
               ENDIF   !EDGE
               ENDIF   !IETOT.GT.1

C             IF (IETOT.LT.NENE.AND.
C     &           EBUFF(IETOT+1).EQ.EBUFF(IETOT)) THEN
C               AMU=ABUFF(IETOT)
C             ENDIF   !EDGE

             COEFTOT(IETOT)=COEFTOT(IETOT)+DENS*THICK*AMU
         ENDDO   !NETOT

      ENDDO !IFILES

      IF (THICKTOT.GT.0.0) THEN
          DENSTOT=DENTHIT/THICKTOT
        IF (DENTHIT.GT.0.0) THEN
          DO IETOT=1,NETOT
             COEFTOT(IETOT)=COEFTOT(IETOT)/DENTHIT
          ENDDO
        ENDIF
      ENDIF

      CLOSE(LUNAM)

      OPEN(UNIT=LUNAMO,FILE=FILEAMO,STATUS='NEW')
        WRITE(LUNAMO,*)CODE
        WRITE(LUNAMO,*)DENSTOT,THICKTOT
        WRITE(LUNAMO,*)NETOT
        DO IETOT=1,NETOT
          WRITE(LUNAMO,*) EBUFFT(IETOT),COEFTOT(IETOT)
        ENDDO
      CLOSE(LUNAMO)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'     Merged absorption coefficients written to file:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     ',FILEAMO
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     Total thickness [m]        :',THICKTOT
      WRITE(LUNGFO,*)'     Effective density [kg/m**3]:',DENSTOT
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

      WRITE(LUNGFO,*)'     Result of merging:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)CODE
      WRITE(LUNGFO,*)DENSTOT,THICKTOT
      WRITE(LUNGFO,*)NETOT
      DO IETOT=1,NETOT
        WRITE(LUNGFO,*) EBUFFT(IETOT),COEFTOT(IETOT)
      ENDDO

      ABSTHI=THICKTOT
      FILEABS=FILEAMO

      RETURN
      END
