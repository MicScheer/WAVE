*CMZ :  4.00/14 01/01/2022  22.11.28  by  Michael Scheer
*CMZ :  4.00/13 20/12/2021  14.08.18  by  Michael Scheer
*CMZ :  4.00/04 17/05/2019  14.22.20  by  Michael Scheer
*CMZ :  3.03/04 02/01/2018  15.10.10  by  Michael Scheer
*CMZ :  3.02/06 15/04/2015  12.10.51  by  Michael Scheer
*CMZ :  3.02/05 22/03/2015  19.10.31  by  Michael Scheer
*CMZ :  3.02/03 10/11/2014  14.50.57  by  Michael Scheer
*CMZ :  3.02/00 24/09/2014  11.50.26  by  Michael Scheer
*CMZ :  3.01/07 23/06/2014  16.07.32  by  Michael Scheer
*CMZ :  3.01/06 20/06/2014  16.51.11  by  Michael Scheer
*CMZ :  3.01/05 13/06/2014  09.06.46  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  10.40.59  by  Michael Scheer
*CMZ :  2.70/11 22/02/2013  14.43.36  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  12.35.44  by  Michael Scheer
*CMZ :  2.68/05 24/09/2012  12.06.38  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  16.02.29  by  Michael Scheer
*CMZ :  2.63/03 10/06/2008  14.12.49  by  Michael Scheer
*CMZ :  2.48/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.48/03 03/03/2004  12.49.39  by  Michael Scheer
*CMZ :  2.42/02 12/09/2002  11.04.58  by  Michael Scheer
*CMZ :  2.41/13 22/08/2002  17.13.40  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.41/05 18/04/2002  11.45.24  by  Michael Scheer
*CMZ :  2.40/02 14/03/2002  16.22.32  by  Michael Scheer
*CMZ :  2.40/00 11/03/2002  16.44.01  by  Michael Scheer
*CMZ :  2.17/00 02/11/2000  16.36.28  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.25.04  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.07.22  by  Michael Scheer
*CMZ :  1.00/00 30/09/97  11.38.28  by  Michael Scheer
*-- Author :    Michael Scheer   24/09/97
      SUBROUTINE MHROUT(JD,ICYCLE,CHOPT)
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

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

      INTEGER NVARMXP,NDIMTP
      PARAMETER (NDIMTP=10)   !SEE ALSO CHARACTER(10) CHTAG
      PARAMETER (NVARMXP=100)

      INTEGER ID,ICYCLE,NID,IDVECT(10000),I,JD,KIND(32),LUNU,ic
      INTEGER IX,IY,JHISASCII,ivar,NVAR,NEVENT,IERR,IJ,II,III,N,NTAGMX
      INTEGER I1DIM,I2DIM,INTUP,ICAL,J,NX,NY,LOC,NWT,INOHEADER
      INTEGER ITIT,JTIT,IGETLASTCHAR
      integer kdcode,istat

      REAL*4 X,Y,XMI,XMA,YMI,YMA,DX,HIm,DY,HIJM
      REAL*4 XNTUP(NVARMXP),RLOW(NVARMXP),RHIGH(NVARMXP)

      CHARACTER(1) CHOPT,C
      CHARACTER(10) CHTAG(NVARMXP),CHDUM
      CHARACTER(16) CRUN,CID
      CHARACTER(80) TIT,TIT1,CHVAR,CLINE
      CHARACTER(128) FILE,FILE1
      character(1024) cwrite
      EXTERNAL IGETLASTCHAR
      equivalence(c,ic)
      data ic/0/

      DATA ICAL/0/
      DATA INOHEADER/0/
      DATA LUNU/99/
      data kdcode/0/

      IF (IHISASCII.NE.0.and.mhbookp.eq.0) THEN

        IF (ICAL.EQ.0) THEN
          JHISASCII=ABS(IHISASCII)
          IF (IHISASCII.LT.0) THEN
            INOHEADER=1
          ENDIF
          I1DIM=(JHISASCII-JHISASCII/10*10)/1
          I2DIM=(JHISASCII-JHISASCII/100*100)/10
          INTUP=(JHISASCII-JHISASCII/1000*1000)/100
          CALL miztoc(ABS(ICODE),CRUN)
          ICAL=1
        ENDIF

        IF (JD.EQ.0) THEN
          CALL hidallm(IDVECT,NID)
        ELSE   !JD.EQ.0
          NID=1
          IDVECT(1)=JD
        ENDIF  !JD.EQ.0

        DO I=1,NID
          ID=IDVECT(I)
          IF (ID.EQ.16) THEN
            LUNGFO=6
          ELSE
            LUNGFO=16
          ENDIF
          call hkindm(ID,KIND,'A')

          CID='          '
          CALL miztoc(ID,CID)
          JTIT=IGETLASTCHAR(1,10,CID,C)

          IF (KIND(1).NE.0) THEN
            IF (I1DIM.NE.0) THEN
              CALL hgivem(ID,TIT,NX,XMI,XMA,NY,YMI,YMA,NWT,LOC)
              DO ITIT=1,IGETLASTCHAR(1,80,TIT,C)
                C=TIT(ITIT:ITIT)
                IF (C.GE.'A'.AND.C.LE.'Z') ic=ic+32
                IF (
     &              C.GE.'A'.AND.C.LE.'Z'
     &              .OR.C.GE.'a'.AND.C.LE.'z'
     &              .OR.C.GE.'0'.AND.C.LE.'9'
     &              .OR.C.EQ.'('.AND.C.EQ.')'
     &              ) THEN
                  TIT1(ITIT:ITIT)=c
                ELSE
                  TIT1(ITIT:ITIT)='_'
                ENDIF
              ENDDO
              ITIT=IGETLASTCHAR(1,80,TIT,C)
              do while (tit1(itit:itit).eq.'_')
                itit=itit-1
              enddo
              FILE=TIT1(1:ITIT)//'_'//CID(1:JTIT)//'.wvh'
              ITIT=IGETLASTCHAR(1,128,FILE,C)
              IF (ITIT.GE.50) THEN
                J=0
                DO II=1,ITIT
                  IF (FILE(II:II).NE.'_') THEN
                    J=J+1
                    FILE1(J:J)=FILE(II:II)
                  ENDIF
                ENDDO
                FILE=FILE1
              ENDIF   !(ITIT.GT.80)

              OPEN(UNIT=LUNU,FILE=FILE)

              IF (INOHEADER.EQ.0) THEN
                WRITE(LUNU,*)CHISASCII//' HBOOK1'
                CLINE(1:1)=CHISASCII
                CLINE(2:2)=' '
                WRITE(CLINE(3:10),'(I8)')ICODE
                CLINE=CLINE(1:10)//' | '//CODE
                WRITE(LUNU,*)CLINE
                WRITE(CLINE(3:10),'(I8)')ID
                CLINE=CLINE(1:10)//' | '//TIT(1:ITIT)
                WRITE(LUNU,*)CLINE
              ENDIF
              DX=(XMA-XMI)/NX
              X=XMI-0.5*DX
              DO IX=1,NX
                X=X+DX
                WRITE(LUNU,*)X,HIm(ID,IX)
              ENDDO   !IX
              CLOSE(LUNU)
            ENDIF !I1DIM
          ELSE IF  (KIND(2).NE.0) THEN
            IF  (I2DIM.NE.0) THEN
              CALL hgivem(ID,TIT,NX,XMI,XMA,NY,YMI,YMA,NWT,LOC)
              DO ITIT=1,IGETLASTCHAR(1,80,TIT,C)
                C=TIT(ITIT:ITIT)
                IF (C.GE.'A'.AND.C.LE.'Z') ic=ic+32
                IF (
     &              C.GE.'A'.AND.C.LE.'Z'
     &              .OR.C.GE.'a'.AND.C.LE.'z'
     &              .OR.C.GE.'0'.AND.C.LE.'9'
     &              .OR.C.EQ.'('.AND.C.EQ.')'
     &              ) THEN
                  TIT1(ITIT:ITIT)=c
                ELSE
                  TIT1(ITIT:ITIT)='_'
                ENDIF
              ENDDO
              ITIT=IGETLASTCHAR(1,80,TIT,C)
              do while (tit1(itit:itit).eq.'_')
                itit=itit-1
              enddo
              FILE=TIT1(1:ITIT)//'_'//CID(1:JTIT)//'.wvh'
              ITIT=IGETLASTCHAR(1,128,FILE,C)
              IF (ITIT.GE.50) THEN
                J=0
                DO II=1,ITIT
                  IF (FILE(II:II).NE.'_') THEN
                    J=J+1
                    FILE1(J:J)=FILE(II:II)
                  ENDIF
                ENDDO
                FILE=FILE1
              ENDIF   !(ITIT.GT.80)

              OPEN(UNIT=LUNU,FILE=FILE)

              IF (INOHEADER.EQ.0) THEN
                WRITE(LUNU,*)CHISASCII//' HBOOK2'
                CLINE(1:1)=CHISASCII
                CLINE(2:2)=' '
                WRITE(CLINE(3:10),'(I8)')ICODE
                CLINE=CLINE(1:10)//' | '//CODE
                WRITE(LUNU,*)CLINE
                WRITE(CLINE(3:10),'(I8)')ID
                CLINE=CLINE(1:10)//' | '//TIT(1:ITIT)
                WRITE(LUNU,*)CLINE
              ENDIF
              DX=(XMA-XMI)/NX
              DY=(YMA-YMI)/NX
              Y=YMI-0.5*DY
              DO IY=1,NY
                Y=Y+DY
                X=XMI-0.5*DX
                DO IX=1,NX
                  X=X+DX
                  WRITE(LUNU,*)X,Y,HIJM(ID,IX,IY)
                ENDDO   !IX
              ENDDO   !IY
              CLOSE(LUNU)
            ENDIF  ! I2DIM
          ELSE IF  (KIND(4).NE.0) THEN
            IF  (INTUP.NE.0) THEN
              IF  (ID.NE.15.AND.ID.NE.16) THEN
                DO IJ=1,NVARMXP
                  CHTAG(IJ)='          '
                ENDDO   !IJ
                CALL HGNPARm(ID,'MHROUT')
                CALL hnoentm(ID,NEVENT)
                NVAR=NVARMXP
                CALL hgivenm(ID,TIT,NVAR,CHTAG,RLOW,RHIGH)
                DO ITIT=1,IGETLASTCHAR(1,80,TIT,C)
                  C=TIT(ITIT:ITIT)
                  IF (C.GE.'A'.AND.C.LE.'Z') ic=ic+32
                  IF (
     &                C.GE.'A'.AND.C.LE.'Z'
     &                .OR.C.GE.'a'.AND.C.LE.'z'
     &                .OR.C.GE.'0'.AND.C.LE.'9'
     &                .OR.C.EQ.'('.AND.C.EQ.')'
     &                ) THEN
                    TIT1(ITIT:ITIT)=c
                  ELSE
                    TIT1(ITIT:ITIT)='_'
                  ENDIF
                ENDDO
                ITIT=IGETLASTCHAR(1,80,TIT,C)
                do while (tit1(itit:itit).eq.'_')
                  itit=itit-1
                enddo
                FILE=TIT1(1:ITIT)//'_'//CID(1:JTIT)//'.wvh'
                ITIT=IGETLASTCHAR(1,128,FILE,C)
                IF (ITIT.GE.50) THEN
                  J=0
                  DO II=1,ITIT
                    IF (FILE(II:II).NE.'_') THEN
                      J=J+1
                      FILE1(J:J)=FILE(II:II)
                    ENDIF
                  ENDDO
                  FILE=FILE1
                ENDIF !(ITIT.GT.80)
                IF (NVAR.GT.NVARMXP) THEN
                  WRITE(6,*)
     &              '*** ERROR IN MHROUT: DIMENSION NVARMXP EXCEEDED ***'
                  WRITE(6,*)'*** NTUPLE: ',ID,TIT
                  STOP
                ENDIF
                NTAGMX=0
                DO 10 II=1,NVARMXP
                  DO IJ=1,NDIMTP
                    CHDUM=CHTAG(II)
                    IF (CHDUM(IJ:IJ).EQ.' ') THEN
                      IF (IJ.GE.NTAGMX) THEN
                        NTAGMX=IJ
                        GOTO 10
                      ENDIF
                    ENDIF
                  ENDDO  !IJ
10              CONTINUE !II

                if (nvar.le.5) then
                  OPEN(UNIT=LUNU,FILE=FILE,recl=256)
                else if (nvar.le.10) then
                  OPEN(UNIT=LUNU,FILE=FILE,recl=512)
                else if (nvar.le.20) then
                  OPEN(UNIT=LUNU,FILE=FILE,recl=1024)
                else if (nvar.le.50) then
                  OPEN(UNIT=LUNU,FILE=FILE,recl=2048)
                else
                  OPEN(UNIT=LUNU,FILE=FILE)
                endif

                IF (INOHEADER.EQ.0) THEN
                  WRITE(LUNU,*)CHISASCII//' hbookn'
                  CLINE(1:1)=CHISASCII
                  CLINE(2:2)=' '
                  WRITE(CLINE(3:10),'(I8)')ICODE
                  CLINE=CLINE(1:10)//' | '//CODE
                  WRITE(LUNU,*)CLINE
                  WRITE(CLINE(3:10),'(I8)')ID
                  CLINE=CLINE(1:10)//' | '//TIT(1:ITIT)
                  WRITE(LUNU,*)CLINE
                  WRITE(LUNU,*)CHISASCII,' ',NVAR,NTAGMX-1,nevent
                  WRITE(LUNU,*)CHISASCII,' ',(CHtag(ivar)(1:ntagmx),ivar=1,nvar)
                  WRITE(LUNU,*)CHISASCII
                ENDIF   !INOHEADER
                DO IX=1,NEVENT
                  CALL hgnfm(ID,IX,XNTUP,IERR)
                  IF (IERR.NE.0) THEN
                    WRITE(LUNGFO,*)
                    WRITE(LUNGFO,*)'*** WARNING IN MHROUT:'
                    WRITE(LUNGFO,*)'ERROR WHILE READING NTUPLE'
                    WRITE(LUNGFO,*)'NTUPLE, EVENT:', ID,IX
                    WRITE(6,*)
                    WRITE(6,*)'*** WARNING IN MHROUT:'
                    WRITE(6,*)'ERROR WHILE READING NTUPLE'
                    WRITE(6,*)'NTUPLE, EVENT:', ID,IX
                  ENDIF  !IERR
                  WRITE(cwrite,*)(XNTUP(N),N=1,NVAR)
                  write(lunu,*)cwrite(2:len_trim(cwrite))
                ENDDO   !NEVENT
                CLOSE(LUNU)
              ELSE IF (INTUP.NE.0.AND.ID.EQ.15) THEN !ID
                OPEN(UNIT=LUNU,FILE=FILE)

                IF (INOHEADER.EQ.0) THEN
                  WRITE(LUNU,*)'! NTUPLE 15 OF WAVE.IN'
                  WRITE(LUNU,*)'!',ICODE,CODE
                  WRITE(LUNU,*)'!'
                ENDIF   !INOHEADER
                REWIND(LUNGFI)
100             READ(LUNGFI,'(A)',END=900)TIT
                WRITE(LUNU,*)TIT
                GOTO 100
900             CLOSE(LUNU)
              ELSE IF (INTUP.NE.0.AND.ID.EQ.16) THEN !ID
                OPEN(UNIT=LUNU,FILE=FILE)
                IF (INOHEADER.EQ.0) THEN
                  WRITE(LUNU,*)'! NTUPLE 16 OF WAVE.OUT'
                  WRITE(LUNU,*)'!',ICODE,CODE
                  WRITE(LUNU,*)'!'
                ENDIF   !INOHEADER
                REWIND(16)
101             READ(16,'(A)',END=901)TIT
                WRITE(LUNU,*)TIT
                GOTO 101
901             CLOSE(LUNU)
              ELSE  !ID
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'*** WARNING IN MHROUT:'
                WRITE(LUNGFO,*)'UNKNOWN HISTOGRAM TYPE'
                WRITE(LUNGFO,*)'CANNOT WRITE ASCII FILE FOR ID',ID
                WRITE(LUNGFO,*)
                WRITE(6,*)
                WRITE(6,*)'*** WARNING IN MHROUT:'
                WRITE(6,*)'UNKNOWN HISTOGRAM TYPE'
                WRITE(6,*)'CANNOT WRITE ASCII FILE FOR ID',ID
                WRITE(6,*)
              ENDIF !ID
            ENDIF !INTUP
          ELSE
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** WARNING IN MHROUT:'
            WRITE(LUNGFO,*)'UNKNOWN HISTOGRAM TYPE'
            WRITE(LUNGFO,*)'CANNOT WRITE ASCII FILE FOR ID',ID
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** WARNING IN MHROUT:'
            WRITE(6,*)'UNKNOWN HISTOGRAM TYPE'
            WRITE(6,*)'CANNOT WRITE ASCII FILE FOR ID',ID
            WRITE(6,*)
          ENDIF   !KIND
        ENDDO  !I
      else IF (IHISASCII.NE.0.and.mhbookp.ne.0) THEN
        call mh_ascii(jd,ihisascii,icode,chisascii,code)
      ENDIF !IHISASCII

      if (iroottrees.ge.0) CALL HROUTm(JD,ICYCLE,CHOPT)

      LUNGFO=16

      if (iroottrees.ne.0) then

      endif !iroottrees

      RETURN
      END
