*CMZ : 00.00/07 21/07/2009  14.26.41  by  Michael Scheer
*CMZ : 00.00/02 26/01/2004  09.37.26  by  Michael Scheer
*-- Author :    Michael Scheer   13/01/2000

        SUBROUTINE UTIL_MAIL(LUNIN,FILEIN,LUNOUT,FILEBAS,ISTATUS)

      IMPLICIT NONE

      INTEGER LUNIN,LUNOUT,ISTATUS,I,I1,I2,MULTI,LENB,IBOUND,LENMS,LENBAS
     &         ,IBLANK

      CHARACTER(1) CQUOTE,CEQUAL,CBLANK
      CHARACTER(80) C80,FILEIN,FILEBAS,BOUNDARY,FILEMS,CEMPTY,FILENAME

      DATA CEQUAL/'='/
      DATA CQUOTE/'"'/
      DATA CBLANK/' '/

      ISTATUS=0
      MULTI=0

      DO I=1,80
          CEMPTY(I:I)=' '
      ENDDO

      I=1
      LENBAS=0
      DO WHILE (FILEBAS(I:I).NE.CBLANK)
          I=I+1
      END DO
      LENBAS=I-1

      OPEN(UNIT=LUNIN,FILE=FILEIN,STATUS='OLD',ERR=99)
      WRITE(6,*)'READING FILE: ',FILEIN(1:80)

1     READ(LUNIN,'(A80)',END=9,ERR=999)C80

      IF (C80(1:35).EQ.'Content-Type: multipart/alternative'
     & .OR. C80(1:29).EQ.'Content-Type: multipart/mixed') THEN

          MULTI=1

10        IBOUND=0
          I1=0
          I2=0

          DO I=1,80
              IF (C80(I:I).EQ.CEQUAL) THEN
                  I1=I+1
              ELSEIF (I.GT.I1.AND.C80(I:I).EQ.CQUOTE) THEN
                  I2=I
                  GOTO 100 !NEXT LINE
              ENDIF
          ENDDO

100       IF (I1.NE.0.AND.I2.NE.0) THEN
             LENB=I2-I1-1
             BOUNDARY=C80(I1+1:I2-1)
          ELSE
                  READ(LUNIN,'(A80)',END=9,ERR=999)C80
                  GOTO 10
          ENDIF

          GOTO 1 !NEXT LINE

      ENDIF   !(C80(1:35).EQ.'Content-Type: multipart/alternative')

      IF (C80(1:24).EQ.'Content-Type: text/plain'
     & .OR. C80(1:24).EQ.'Content-Type: TEXT/PLAIN') THEN

11        READ(LUNIN,'(A80)',END=9,ERR=999)C80

          FILENAME=FILEBAS(1:LENBAS)//'.TXT'

          IF (C80(1:6).EQ.' name=') THEN
          FILENAME=CEMPTY
          DO I=8,80
              IF (C80(I:I).EQ.CQUOTE) THEN
                  GOTO 1119
              ELSE
                  FILENAME(I-7:I-7)=C80(I:I)
              ENDIF
          ENDDO
          ENDIF

C         IF (C80(1:26).NE.'Content-Transfer-Encoding:') GOTO 11
          IF (C80(1:26).NE.CEMPTY) GOTO 11

1119      OPEN(UNIT=LUNOUT,FILE=FILENAME,STATUS='NEW'
     &      )
          WRITE(6,*)'WRITING FILE: ',FILENAME

111       READ(LUNIN,'(A80)',END=9,ERR=999)C80

          IF (LENB.GT.0) THEN

             IBOUND=0
             DO I=1,80-LENB
                 IF (C80(I:LENB+I-1).EQ.BOUNDARY) THEN
                     IBOUND=1
                     GOTO 1111
                 ENDIF
             ENDDO

          ENDIF   !(LENB.GT.0)

1111      IF (IBOUND.EQ.0) THEN
              WRITE(LUNOUT,'(A80)')C80
              GOTO 111
          ELSE
              CLOSE(LUNOUT)
          ENDIF

          GOTO 1  !NEXT LINE

      ENDIF   !text/plain

      IF (C80(1:23).EQ.'Content-Type: text/html'
     & .OR. C80(1:23).EQ.'Content-Type: TEXT/HTML') THEN

22        READ(LUNIN,'(A80)',END=9,ERR=999)C80

          IF (C80(1:26).NE.'Content-Transfer-Encoding:') GOTO 22

          OPEN(UNIT=LUNOUT,FILE=FILEBAS(1:LENBAS)//'.HTML',STATUS='NEW'
     &      )
          WRITE(6,*)'WRITING FILE: ',FILEBAS(1:LENBAS)//'.HTML'

222       READ(LUNIN,'(A80)',END=9,ERR=999)C80

          IBOUND=0
          DO I=1,80-LENB
              IF (C80(I:LENB+I-1).EQ.BOUNDARY) THEN
                  IBOUND=1
                  GOTO 2222
              ENDIF
          ENDDO

2222      IF (IBOUND.EQ.0) THEN
              WRITE(LUNOUT,'(A80)')C80
              GOTO 222
          ELSE
              CLOSE(LUNOUT)
          ENDIF

          GOTO 1  !NEXT LINE

      ENDIF   !text/html

      IF (C80(1:32).EQ.'Content-Type: application/msword'
     & .OR. C80(1:29).EQ.'Content-Type: application/rtf'
     & .OR. C80(1:29).EQ.'Content-Type: application/pdf'
     & ) THEN

331       I1=0
          I2=0

          DO I=1,80
              IF (I1.EQ.0.AND.C80(I:I).EQ.CQUOTE) THEN
                  I1=I
              ELSEIF (C80(I:I).EQ.CQUOTE) THEN
                  I2=I
              ENDIF
          ENDDO

          IF (I1.EQ.0.OR.I2.EQ.0) THEN
              READ(LUNIN,'(A80)',END=9,ERR=999)C80
              GOTO 331 !NEXT LINE
          ENDIF

          LENMS=I2-1-I1
          FILEMS=C80(I1+1:I2-1)

33        READ(LUNIN,'(A80)',END=9,ERR=999)C80

          IF (C80(1:26).NE.'Content-Transfer-Encoding:') GOTO 33

            IF (C80(28:33).NE.'base64') THEN
              WRITE(6,*)
              WRITE(6,*)'*** WARNING IN UTIL_MAIL: UNKNOWN ENCODING-TYPE'
              WRITE(6,*)C80
              WRITE(6,*)'    DETECTED ***'
              WRITE(6,*)
              ISTATUS=ISTATUS+1
              GOTO 1
          ENDIF   !BASE64

          OPEN(UNIT=LUNOUT,FILE=FILEMS(1:LENMS)//'_B64',STATUS='NEW'
     &      )
          WRITE(6,*)'WRITING FILE: ',FILEMS(1:LENMS)//'_B64'

C LOOK FOR BLANK LINE

332       READ(LUNIN,'(A80)',END=9,ERR=999)C80

          IBLANK=1
          DO I=1,80
              IF (C80(I:I).NE.CBLANK) IBLANK=0
          ENDDO

          IF (IBLANK.EQ.0) GOTO 332

333       READ(LUNIN,'(A80)',END=9,ERR=999)C80

          IBOUND=0
          DO I=1,80-LENB
              IF (C80(I:LENB+I-1).EQ.BOUNDARY) THEN
                  IBOUND=1
                  GOTO 3333
              ENDIF
          ENDDO

3333      IF (IBOUND.EQ.0) THEN
              DO I=1,80
                  IF (C80(I:I).EQ.CBLANK) GOTO 123
              ENDDO
123           IF (I.GT.1) THEN
                  WRITE(LUNOUT,'(A)')C80(1:I-1)
              ENDIF
              GOTO 333
          ELSE
              CLOSE(LUNOUT)
          ENDIF

          GOTO 1  !NEXT LINE

      ENDIF   !application/msword

      IF (C80(1:36).EQ.'Content-Type: application/postscript'
     & ) THEN

441       I1=0
          I2=0

          DO I=1,80
              IF (I1.EQ.0.AND.C80(I:I).EQ.CQUOTE) THEN
                  I1=I
              ELSEIF (C80(I:I).EQ.CQUOTE) THEN
                  I2=I
              ENDIF
          ENDDO

          IF (I1.EQ.0.OR.I2.EQ.0) THEN
              READ(LUNIN,'(A80)',END=9,ERR=999)C80
              GOTO 441 !NEXT LINE
          ENDIF

          LENMS=I2-1-I1
          FILEMS=C80(I1+1:I2-1)

44        READ(LUNIN,'(A80)',END=9,ERR=999)C80

          IF (C80(1:26).NE.'Content-Transfer-Encoding:') GOTO 44

          OPEN(UNIT=LUNOUT,FILE=FILEMS(1:LENMS),STATUS='NEW'
     &      )
          WRITE(6,*)'WRITING FILE: ',FILEMS(1:LENMS)

C LOOK FOR BLANK LINE

442       READ(LUNIN,'(A80)',END=9,ERR=999)C80

          IBLANK=1
          DO I=1,80
              IF (C80(I:I).NE.CBLANK) IBLANK=0
          ENDDO

          IF (IBLANK.EQ.0) GOTO 442

444       READ(LUNIN,'(A80)',END=9,ERR=999)C80

          IF (MULTI.NE.0) THEN

             IBOUND=0
             DO I=1,80-LENB
                 IF (C80(I:LENB+I-1).EQ.BOUNDARY) THEN
                     IBOUND=1
                     GOTO 4444
                 ENDIF
             ENDDO

          ENDIF   !MULTI

4444      IF (IBOUND.EQ.0) THEN
              DO I=80,1,-1
                  IF (C80(I:I).NE.CBLANK) GOTO 124
              ENDDO
124           IF (I.GT.1) THEN
                  WRITE(LUNOUT,'(A)')C80(1:I)
              ENDIF
              GOTO 444
          ELSE
              CLOSE(LUNOUT)
          ENDIF

          GOTO 1  !NEXT LINE

      ENDIF   !application/postscript

      IF (C80(1:13).EQ.'Content-Type:') THEN
          WRITE(6,*)
          WRITE(6,*)'*** WARNING IN UTIL_MAIL: UNKNOWN CONTENT-TYPE'
          WRITE(6,*)C80
          WRITE(6,*)'    DETECTED ***'
          WRITE(6,*)
          ISTATUS=ISTATUS+1
      ENDIF

      GOTO 1  !NEXT LINE

9     CLOSE(LUNIN)
      RETURN

99    ISTATUS=-1
      WRITE(6,*)'*** ERROR DURING OPENING IN UTIL_MAIL ***'
      RETURN

999   ISTATUS=-2
      WRITE(6,*)'*** ERROR DURING READ IN UTIL_MAIL ***'
      RETURN

      END
